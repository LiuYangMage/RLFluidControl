import tensorflow as tf
import numpy as np
import random
import os
from utils import ReplayBuffer

os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"  
os.environ["CUDA_VISIBLE_DEVICES"]= str(1)

class Agent:
    def __init__(self, state_dim, action_dim, explore_noise = "Gaussian", *args, **kwargs):
        self.lr = 1e-4
        self.gamma = 0.99
        self.tau = 0.005
        self.bs = 512
        self.bfs = 1000000
        self.d = 2
        self.explore_noise = explore_noise
        self.explore_noise_size = 0.1
        self.criticreg_noise_size = 0.2
        self.criticreg_noise_clip = 0.5

        self.state_dim = state_dim
        self.action_dim = action_dim
        self.actor_nn_dim = [256, 256, self.action_dim]
        self.critic_nn_dim = [256, 256, 1]

        self.state1_place = tf.placeholder(tf.float32, [None, self.state_dim])
        self.action_place = tf.placeholder(tf.float32, [None, self.action_dim])
        self.reward_place = tf.placeholder(tf.float32, [None,1])
        self.isdone_place = tf.placeholder(tf.float32, [None,1])
        self.state2_place = tf.placeholder(tf.float32, [None, self.state_dim])


        with tf.variable_scope("target_actor", reuse = tf.AUTO_REUSE):
            self.Q_next_action = self.actor_nn(self.state2_place)
        self.Q_next_noise = tf.clip_by_value(tf.random.normal([self.bs, self.action_dim], 0, self.criticreg_noise_size),
                         - self.criticreg_noise_clip, self.criticreg_noise_clip)
        self.Q_next_noisy_action = tf.clip_by_value(self.Q_next_action + self.Q_next_noise, -1, 1)
        with tf.variable_scope("target_critic_1", reuse = tf.AUTO_REUSE):
            self.Q_critic_1 = self.critic_nn(self.state2_place, self.Q_next_noisy_action)
        with tf.variable_scope("target_critic_2", reuse = tf.AUTO_REUSE):
            self.Q_critic_2 = self.critic_nn(self.state2_place, self.Q_next_noisy_action)
        self.Q_critic_min = tf.minimum(self.Q_critic_1, self.Q_critic_2)
        self.Q_y = self.reward_place + self.gamma * (1-self.isdone_place) * self.Q_critic_min
        with tf.variable_scope("main_critic_1", reuse = tf.AUTO_REUSE):
            self.Q_Q_1 = self.critic_nn(self.state1_place, self.action_place)
        with tf.variable_scope("main_critic_2", reuse = tf.AUTO_REUSE):
            self.Q_Q_2 = self.critic_nn(self.state1_place, self.action_place) 
        self.Q_loss = tf.reduce_mean((self.Q_Q_1 - self.Q_y)**2) + tf.reduce_mean((self.Q_Q_2 - self.Q_y)**2)

        with tf.variable_scope("main_actor", reuse = tf.AUTO_REUSE):
            self.P_this_action = self.actor_nn(self.state1_place)
        with tf.variable_scope("main_critic_1", reuse = tf.AUTO_REUSE):
            self.P_Q_1 = self.critic_nn(self.state1_place, self.P_this_action)
        self.P_loss = - tf.reduce_mean(self.P_Q_1)

        with tf.variable_scope("main_actor", reuse = tf.AUTO_REUSE):
            self.action = self.actor_nn(self.state1_place)


        all_variables = tf.trainable_variables()
        self.main_critic_var = [i for i in all_variables if "main_critic" in i.name]
        self.target_critic_var = [i for i in all_variables if "target_critic" in i.name]
        self.main_actor_var = [i for i in all_variables if "main_actor" in i.name]
        self.target_actor_var = [i for i in all_variables if "target_actor" in i.name]

        assert len(self.main_critic_var) == len(self.target_critic_var)
        assert len(self.main_actor_var) == len(self.target_actor_var)

        self.Q_op =  tf.train.AdamOptimizer(self.lr).minimize(self.Q_loss, var_list = self.main_critic_var) 
        self.P_op =  tf.train.AdamOptimizer(self.lr).minimize(self.P_loss, var_list = self.main_actor_var) 
        self.T_init = [tf.assign(T, M) for (T,M) in zip(self.target_critic_var + self.target_actor_var, 
                                                        self.main_critic_var + self.main_actor_var)]
        self.T_op = [tf.assign(T, self.tau * M + (1 - self.tau) * T) for (T,M) in zip(
            self.target_critic_var + self.target_actor_var, self.main_critic_var + self.main_actor_var)]


        self.replay_buffer = ReplayBuffer(self.state_dim, self.action_dim, self.bfs)

        self.step_count = 0
        self.total_step_count = 0
        self.episode_count = 0
        self.train_count = 0

        config = tf.ConfigProto()
        config.gpu_options.allow_growth = True
        self.sess = tf.Session(config=config)
        self.sess.run(tf.global_variables_initializer())
        self.sess.run(self.T_init)
        self.saver = tf.train.Saver(max_to_keep=1000)


    def actor_nn(self, state, bound = True):
        dim = self.actor_nn_dim
        A = state
        for i in range(0,len(dim)-1):
            A = tf.layers.dense(A, units= dim[i], activation = tf.nn.relu)
        action = tf.layers.dense(A, units= dim[-1], activation = tf.nn.tanh)
        return action


    def critic_nn(self, state, action):
        dim = self.critic_nn_dim
        A = tf.concat([state, action], axis = 1)
        for i in range(0,len(dim)-1):
            A = tf.layers.dense(A, units= dim[i], activation = tf.nn.relu)
        critic = tf.layers.dense(A, units= dim[-1], activation = None)
        return critic



    def get_action(self, state_data, stochastic = True):
        this_action = self.sess.run(self.action, feed_dict= {self.state1_place: state_data})
        if stochastic:
            if self.explore_noise == "Gaussian":
                explore_noise = np.random.normal(0, self.explore_noise_size, [1, self.action_dim])
            else:
                raise NotImplementedError
            this_action = np.clip(this_action + explore_noise, -1, 1)
        return this_action

    def eval_loss(self, bs = None):
        if bs is None:
            bs = self.bs
        this_bs = np.minimum(bs, self.replay_buffer.size)
        this_batch = self.replay_buffer.sample_batch(this_bs)
        feed_dict = {self.state1_place: this_batch["obs1"],
                     self.action_place: this_batch["acts"],
                     self.reward_place: this_batch["rews"],
                     self.isdone_place: this_batch["done"],
                     self.state2_place: this_batch["obs2"]}
        pass

    def train_iter(self):

        if self.bs <= self.replay_buffer.size:
            this_batch = self.replay_buffer.sample_batch(self.bs)
            
            feed_dict = {self.state1_place: this_batch["obs1"],
                        self.action_place: this_batch["acts"],
                        self.reward_place: this_batch["rews"],
                        self.isdone_place: this_batch["done"],
                        self.state2_place: this_batch["obs2"]}

            self.sess.run([self.Q_op], feed_dict=feed_dict)
            if self.total_step_count % self.d == 0:
                self.sess.run([self.P_op], feed_dict=feed_dict)
                self.sess.run(self.T_op)

            self.train_count += 1
            self.total_step_count += 1


    def record(self, this_state, this_action, this_reward, this_done, next_state):
        self.replay_buffer.store(obs=this_state, 
                                act=this_action,
                                rew=this_reward,
                                next_obs=next_state,
                                done=this_done)
        self.step_count += 1


    def reset_agent(self):

        self.replay_buffer = ReplayBuffer(self.state_dim, self.action_dim, self.bfs)

        self.step_count = 0
        self.total_step_count = 0
        self.train_count = 0
        self.episode_count = 0

        self.sess.run(tf.global_variables_initializer())
        self.sess.run(self.T_init)

    def reset_episode(self):

        self.step_count = 0
        self.train_count = 0
        self.episode_count += 1

