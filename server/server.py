import os
import pickle
from xmlrpc.server import SimpleXMLRPCServer
import tensorflow as tf
import numpy as np
from agent_TD3 import Agent
from utils import KalmanFilter, NoneFilter
from datetime import datetime
import argparse




class Server():
    def __init__(self, environment, filter, explore_noise, *args, **kwargs):
        self.environment = environment
        self.agent = Agent(2, 2, explore_noise = explore_noise)
        self._stamp("Environment: " + environment)
        self._stamp("Using Exploration Noise: " + explore_noise)
        self.state_record = [] # filtered.
        self.unfiltered_state_record = []
        self.action_record = []
        self.state_mean = np.zeros((1,2))
        # self.server = SimpleXMLRPCServer(("0.0.0.0", 8000))
        self.server = SimpleXMLRPCServer(("localhost", 8000),logRequests=False)
        self._register()
        self.filter_name = filter
        if self.filter_name == "Kalman":
            self.filter = KalmanFilter()
            self._stamp("Using Kalman Filter")
        elif self.filter_name == "None":
            self.filter = NoneFilter()
            self._stamp("Using No Filter")
        else:
            raise NotImplementedError
        self.save_model_dir = "save"
        self.save_data_dir = "save_data"
        self.save_eval_dir = "save_eval"


    def _stamp(self, string):
        time = "UTC " + datetime.utcnow().isoformat(sep=" ", timespec="milliseconds") + " "
        print(time + string, flush = True)

    def _register(self):
        self.server.register_function(self._init, "init")
        self.server.register_function(self._start_episode, "start_episode")
        self.server.register_function(self._request_stochastic_action, "request_stochastic_action")
        self.server.register_function(self._request_deterministic_action, "request_deterministic_action")
        self.server.register_function(self._train, "train")
        self.server.register_function(self._save, "save")
        self.server.register_function(self._save_eval, "save_eval")
        self.server.register_function(self._restore, "restore")

    def _reward_func(self, this_state, action, next_state):
        reward = -(np.sign(next_state[0,1]) * (next_state[0,1])**2)
        #reward = -(np.sign(next_state[0,1]) * ((next_state[0,1])**2) + 0.1*((next_state[0,0])**2))
        return reward

    def _converter(self, signal):
        """
        convert the signal from sensor to forces and moment.
        Args: 
            signal list: 6*n, 
        Return:
            states (ndarray): (n,2), dimensionless C_lift and C_drag 
        """
        non_dimensional = 0.5 * 1000 * 0.2**2 * 2 * 20 * 0.0254**2
        calmat = np.array([
        [-0.03424, 0.00436,  0.04720,  -1.87392,  0.00438,  1.89084],
        [-0.02353,  2.19709,  -0.02607,  -1.08484,  0.03776,  -1.09707],
        [3.32347,  -0.13564,  3.36623,  -0.09593,  3.35105,  -0.15441],
        [-0.00778,  1.04073,  -3.83509,  -0.40776,  3.81653,  -0.69079],
        [4.36009,  -0.17486,  -2.26568,  0.94542,  -2.18742,  -0.79758],
        [0.03483,  -2.32692,  0.02473,  -2.30260,  0.01156,  -2.33458]])

        signal_array = np.mean(np.array(signal).reshape((6,-1)), axis = 1, keepdims = True)
        calforce = calmat @ signal_array
        lift = calforce[0:1,:]*4.44822 / non_dimensional # C_l
        drag = calforce[1:2,:]*4.44822 / non_dimensional # C_d
        state = np.concatenate([lift, drag], axis = 0) #(2, 1)

        return state.transpose() #(1,2)

    def _init(self, episode_count):
        try:
            self._restore(episode_count)
            return True

        except Exception:
            self.agent.reset_agent()
            self._stamp("Initialized!")
            return False

    def _start_episode(self,raw_data):
        # raw_data for calibrating
        if self.environment == "Experiment":
            self.state_mean = self._converter(raw_data) #(1, 2)
        self.state_record = []
        self.unfiltered_state_record = []
        self.action_record = []
        self.agent.reset_episode()
        self.filter.reset()
        self._stamp("Epsode Start!")
        return True

    def _request_action(self, raw_data, stochastic):
        if self.environment == "Experiment":
            calibrated_state = self._converter(raw_data) - self.state_mean #(1,2)
        elif self.environment == "CFD":
            calibrated_state = np.asfarray(raw_data.split("_"),float)[None,:]
        else:
            raise NotImplementedError
        filtered_state = self.filter.estimate(calibrated_state)

        action = self.agent.get_action(filtered_state, stochastic = stochastic) #(1,2)
        self.unfiltered_state_record.append(calibrated_state)
        self.state_record.append(filtered_state)
        self.action_record.append(action)
        self._stamp("C_lift: {:.3f} C_drag: {:.3f} Action: {:.3f} {:.3f} ".format(
                    filtered_state[0,0], 
                    filtered_state[0,1],
                    action[0,0],
                    action[0,1]))

        if self.environment == "Experiment":
            raw_action = action[0,:].tolist()
        elif self.environment == "CFD":
            raw_action = str(action[0,0])+"_"+str(action[0,1])
        else:
            raise NotImplementedError
        return raw_action

    def _request_stochastic_action(self, raw_data):
        return self._request_action(raw_data, True)

    def _request_deterministic_action(self, raw_data):
        return self._request_action(raw_data, False)
       
    def _train(self, steps):
        if not os.path.exists(self.save_data_dir):
            os.mkdir(self.save_data_dir)
        np.savez(self.save_data_dir + "/data_{}.npz".format(self.agent.episode_count), 
                    state = np.array(self.state_record), 
                    unfiltered_state = np.array(self.unfiltered_state_record), 
                    action = np.array(self.action_record))
        record_length = len(self.state_record)
        self._stamp("Length of Record: " + str(record_length))
        for i in range(record_length - 1):
            reward = self._reward_func(self.state_record[i], 
                                       self.action_record[i], 
                                       self.state_record[i+1])
            self.agent.replay_buffer.store(self.state_record[i],
                                        self.action_record[i], 
                                        reward,
                                        self.state_record[i+1],
                                        0)
        self._stamp("Training Start!")        
        for i in range(steps):
            self.agent.train_iter()
        self._stamp("Training End!")        
        return True

    def _save_eval(self, episode_count, rep_count):
        if not os.path.exists(self.save_eval_dir):
            os.mkdir(self.save_eval_dir)
        np.savez(self.save_eval_dir + "/data_{}_{}.npz".format(episode_count, rep_count), 
            state = np.array(self.state_record), 
            unfiltered_state = np.array(self.unfiltered_state_record), 
            action = np.array(self.action_record))
        return True

    def _save(self, dummy = None):
        if not os.path.exists(self.save_model_dir):
            os.mkdir(self.save_model_dir)
        self.agent.saver.save(self.agent.sess, self.save_model_dir + "/{}.ckpt".format(self.agent.episode_count))
        
        pickle_out = open(self.save_model_dir + "/{}.pickle".format(self.agent.episode_count),"wb")
        pickle.dump(self.agent.replay_buffer, pickle_out)
        pickle_out.close()
        self._stamp("Saved Episode {}!".format(self.agent.episode_count))
        return True

    def _restore(self,episode_count):
        self.agent.saver.restore(self.agent.sess, self.save_model_dir + "/{}.ckpt".format(episode_count))
        pickle_in = open(self.save_model_dir + "/{}.pickle".format(episode_count),"rb")
        self.agent.replay_buffer = pickle.load(pickle_in)
        self.agent.episode_count = episode_count
        self._stamp("Restored from Episode {}!".format(episode_count))
        return True

    def start_server(self):
        self._stamp("Server Listening...")
        self.server.serve_forever()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="server for experiments and CFD")
    parser.add_argument("-env", "--environment", choices=["CFD", "Experiment"], default="Experiment")
    parser.add_argument("-fil", "--filter", choices=["None", "Kalman"], default="Kalman")
    parser.add_argument("-exn", "--explore_noise", choices=["Gaussian"], default="Gaussian")
    args = parser.parse_args()

    myserver = Server(**args.__dict__)
    myserver.start_server()