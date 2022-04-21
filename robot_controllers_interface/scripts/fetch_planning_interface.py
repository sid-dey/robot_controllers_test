#!/usr/bin/env python

import sys
import copy
import actionlib
import rospy
import numpy as np

from math import sin, cos

from geometry_msgs.msg import TwistStamped, PoseStamped
from std_msgs.msg import Int8, Float32
from sensor_msgs.msg import JointState
from control_msgs.msg import FollowJointTrajectoryAction, FollowJointTrajectoryGoal
from control_msgs.msg import GripperCommandAction, GripperCommandGoal
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint

DES_JOINTS_TOPIC = "/fetch/des_states"

NO_DOF = 7

class FollowTrajectoryClient(object):

    def __init__(self, name, joint_names):
        self.trajectory = JointTrajectory()
        self.sub_des_states = rospy.Subscriber(DES_JOINTS_TOPIC, JointTrajectory, self.handle_des_states, tcp_nodelay=True,
                                           queue_size=1)    

        self.client = actionlib.SimpleActionClient("%s" % name,
            FollowJointTrajectoryAction)    
        
        rospy.loginfo("Waiting for %s..." % name)
        self.client.wait_for_server()        

        self.trajectory = JointTrajectory()
        self.trajectory.joint_names = joint_names

        point_msg = JointTrajectoryPoint()
        point_msg.positions = np.zeros(NO_DOF) # np.zeros(NO_DOF)
        point_msg.velocities = np.zeros(NO_DOF)
        point_msg.time_from_start = rospy.Duration(2.0)
        self.trajectory.points = [point_msg]

        self.move_to()

        # set up flag for new trajectory received
        self.new_trajectory_received_flag = False

    def handle_des_states(self, data):
        self.trajectory.points = data.points
        self.new_trajectory_received_flag = True

    def move_to(self):
        follow_goal = FollowJointTrajectoryGoal()
        follow_goal.trajectory = self.trajectory
        self.new_trajectory_received_flag = False
        
        self.client.send_goal(follow_goal)
        # self.client.wait_for_result()

if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("usage: fetch_planning_interface.py <planner_name>")
        exit(-1)
    
    # Create a node
    rospy.init_node("fetch_planning_interface")

    # Make sure sim time is working
    while not rospy.Time.now():
        pass

    # Setup clients

    arm_action = FollowTrajectoryClient(sys.argv[1], ["shoulder_pan_joint", "shoulder_lift_joint", \
                                                             "upperarm_roll_joint", "elbow_flex_joint", \
                                                             "forearm_roll_joint", "wrist_flex_joint", \
                                                             "wrist_roll_joint"])

    rospy.loginfo("Moving the arm with position commands")

    while not rospy.is_shutdown():
        if arm_action.new_trajectory_received_flag:
            arm_action.move_to()

