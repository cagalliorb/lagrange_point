#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from math import *
import numpy as np

G = 6.672 * pow(10, -11)  # 万有引力定数
ME = 5.972 * pow(10, 24)  # 地球の質量[kg]
MM = 7.346 * pow(10, 22)  # 月の質量[kg]
MS = 10  # 衛星の質量[kg]
M = ME + MM  # 地球と月の質量の和
R = 3.844 * pow(10, 8)  # 地球と月の距離[m]
OMEGA = sqrt(G * M / (R ** 3)) / (2 * pi)  # 地球と月の重心周りの角速度(ただし単位を[rad/s]->[1/s]とする)
L1 = R * MM / M  # 地球と重心の距離[m]
L2 = R * ME / M  # 月と重心の距離[m]
A = (R ** 3) * (OMEGA ** 2) / (G * MM)  # 無次元化パラメータ1
B = (R ** 3) * (OMEGA ** 2) / (G * ME)  # 無次元化パラメータ2
C = (R ** 3) * (OMEGA ** 2) / (G * MS)  # 無次元化パラメータ3


def rot_mat(theta, all_list):
    tlist = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    return np.dot(tlist, all_list)


def get_dvxe_dt(all_list):
    return 1 / A * (all_list[6] - all_list[0]) / ((all_list[6] - all_list[0]) ** 2 + (
        all_list[7] - all_list[1]) ** 2 + (all_list[8] - all_list[2]) ** 2) ** (3 / 2) \
           + 1 / C * (all_list[12] - all_list[0]) / ((all_list[12] - all_list[0]) ** 2 + (all_list[13] - all_list[1])
                                                     ** 2 + (all_list[14] - all_list[2]) ** 2) ** (3 / 2)


def get_dvye_dt(all_list):
    return 1 / A * (all_list[7] - all_list[1]) / (
                                                     (all_list[6] - all_list[0]) ** 2 + (
                                                         all_list[7] - all_list[1]) ** 2 + (
                                                         all_list[8] - all_list[2]) ** 2) ** (3 / 2) \
           + 1 / C * (all_list[13] - all_list[1]) / ((all_list[12] - all_list[0]) ** 2 + (all_list[13] - all_list[1])
                                                     ** 2 + (all_list[14] - all_list[2]) ** 2) ** (3 / 2)


def get_dvze_dt(all_list):
    return 1 / A * (all_list[8] - all_list[2]) / (
                                                     (all_list[6] - all_list[0]) ** 2 + (
                                                         all_list[7] - all_list[1]) ** 2 + (
                                                         all_list[8] - all_list[2]) ** 2) ** (3 / 2) \
           + 1 / C * (all_list[14] - all_list[2]) / ((all_list[12] - all_list[0]) ** 2 + (all_list[13] - all_list[1])
                                                     ** 2 + (all_list[14] - all_list[2]) ** 2) ** (3 / 2)


def get_dvxm_dt(all_list):
    return 1 / B * (all_list[0] - all_list[6]) / (
                                                     (all_list[6] - all_list[0]) ** 2 + (
                                                         all_list[7] - all_list[1]) ** 2 + (
                                                         all_list[8] - all_list[2]) ** 2) ** (3 / 2) \
           + 1 / C * (all_list[12] - all_list[6]) / ((all_list[6] - all_list[12]) ** 2 + (all_list[7] - all_list[13])
                                                     ** 2 + (all_list[8] - all_list[14]) ** 2) ** (3 / 2)


def get_dvym_dt(all_list):
    return 1 / B * (all_list[1] - all_list[7]) / (
                                                     (all_list[6] - all_list[0]) ** 2 + (
                                                         all_list[7] - all_list[1]) ** 2 + (
                                                         all_list[8] - all_list[2]) ** 2) ** (3 / 2) \
           + 1 / C * (all_list[13] - all_list[7]) / ((all_list[6] - all_list[12]) ** 2 + (all_list[7] - all_list[13])
                                                     ** 2 + (all_list[8] - all_list[14]) ** 2) ** (3 / 2)


def get_dvzm_dt(all_list):
    return 1 / B * (all_list[2] - all_list[8]) / (
                                                     (all_list[6] - all_list[0]) ** 2 + (
                                                         all_list[7] - all_list[1]) ** 2 + (
                                                         all_list[8] - all_list[2]) ** 2) ** (3 / 2) \
           + 1 / A * (all_list[14] - all_list[8]) / ((all_list[12] - all_list[6]) ** 2 + (all_list[13] - all_list[7])
                                                     ** 2 + (all_list[14] - all_list[8]) ** 2) ** (3 / 2)


def get_dvxs_dt(all_list):
    return 1 / B * (all_list[0] - all_list[12]) / ((all_list[12] - all_list[0]) ** 2 + (
        all_list[13] - all_list[1]) ** 2 + (
                                                       all_list[14] - all_list[2]) ** 2) ** (3 / 2) \
           + 1 / A * (all_list[6] - all_list[12]) / ((all_list[12] - all_list[6]) ** 2 + (all_list[13] - all_list[7])
                                                     ** 2 + (all_list[14] - all_list[8]) ** 2) ** (3 / 2)


def get_dvys_dt(all_list):
    return 1 / B * (all_list[1] - all_list[13]) / ((all_list[12] - all_list[0]) ** 2 + (
        all_list[13] - all_list[1]) ** 2 + (
                                                       all_list[14] - all_list[2]) ** 2) ** (3 / 2) \
           + 1 / A * (all_list[7] - all_list[13]) / ((all_list[12] - all_list[6]) ** 2 + (all_list[13] - all_list[7])
                                                     ** 2 + (all_list[14] - all_list[8]) ** 2) ** (3 / 2)


def get_dvzs_dt(all_list):
    return 1 / B * (all_list[2] - all_list[14]) / ((all_list[12] - all_list[0]) ** 2 + (
        all_list[13] - all_list[1]) ** 2 + (
                                                       all_list[14] - all_list[2]) ** 2) ** (3 / 2) \
           + 1 / A * (all_list[8] - all_list[14]) / ((all_list[12] - all_list[6]) ** 2 + (all_list[13] - all_list[7])
                                                     ** 2 + (all_list[14] - all_list[8]) ** 2) ** (3 / 2)


def make_derivative(all_list):
    return np.array([all_list[3], all_list[4], all_list[5],
                     get_dvxe_dt(all_list), get_dvye_dt(all_list), get_dvze_dt(all_list),
                     all_list[9], all_list[10], all_list[11],
                     get_dvxm_dt(all_list), get_dvym_dt(all_list), get_dvzm_dt(all_list),
                     all_list[15], all_list[16], all_list[17],
                     get_dvxs_dt(all_list), get_dvys_dt(all_list), get_dvzs_dt(all_list)])


def cal_and_plot(tend):
    dlt = 0.001
    tlist = np.arange(0.0, tend, dlt)
    lis_num = tlist.shape[0]
    # [xe, ye, ze, vxe, vye, vze, xm, ym, zm, vxm, vym, vzm, xs, ys, zs, vxs, vys, vzs]のn*18配列を用意
    all_list = np.zeros((lis_num, 18))
    xs = np.zeros((lis_num, 3))
    xs_rel = np.zeros((lis_num, 3))
    # 初期条件を代入
    all_list[0, :] = [-L1 / R, 0, 0,
                      0, -L1 / R * (2 * pi), 0,
                      L2 / R, 0, 0,
                      0, L2 / R * (2 * pi), 0,
                      (L2 - L1) / (2 * R), sqrt(3) / 2, 0,
                      -sqrt(3) / 2 * (2 * pi), 1 / 2 * (2 * pi), 0]
    lag_4 = np.array([[(L2 - L1) / (2 * R), sqrt(3) / 2, 0] for n in range(lis_num)])

    xs[0, :] = all_list[0, 12: 15]
    theta = atan2(all_list[0][7], all_list[0][6])
    xs[0, :] = rot_mat(-theta, xs[0][:])
    xs_rel[0][:] = xs[0][0:3] - lag_4[0][0:3]  # L4からの相対座標の計算

    """RK4による計算"""
    for n in range(lis_num - 1):
        k1 = make_derivative(all_list[n])
        k2 = make_derivative(all_list[n] + dlt / 2 * k1)
        k3 = make_derivative(all_list[n] + dlt / 2 * k2)
        k4 = make_derivative(all_list[n] + dlt * k3)
        all_list[n + 1, :] = all_list[n, :] + dlt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        xs[n + 1, :] = all_list[n + 1, 12: 15]
        theta = atan2(all_list[n + 1][7], all_list[n + 1][6])
        xs[n + 1][:] = rot_mat(-theta, xs[n + 1][:])
        xs_rel[n + 1][:] = xs[n + 1][0: 3] - lag_4[n + 1][0: 3]
        print("division: {0}/{1}".format(n, lis_num))

    fig = plt.figure()
    ax1 = fig.add_subplot(111, aspect=1)  # 軌道のxy平面投影グラフの作成
    ax1.plot(all_list[:, 0], all_list[:, 1], label="earth", color="blue")
    ax1.plot(all_list[:, 6], all_list[:, 7], label="moon", color="black")
    ax1.plot(all_list[:, 12], all_list[:, 13], label="satellite", color="red")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_xlim([-1.5, 1.5])
    ax1.set_ylim([-1.5, 1.5])
    plt.legend()
    plt.show()

    """L4に対する相対位置のplot"""
    """
    fig = plt.figure()

    ax1 = fig.add_subplot(131, aspect=1)  # 軌道のxy平面投影グラフの作成
    ax1.plot(xs_rel[:, 0], xs_rel[:, 1])
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_xlim([-0.2, 0.2])
    ax1.set_ylim([-0.2, 0.2])

    ax2 = fig.add_subplot(132, aspect=1)  # 軌道のyz平面投影グラフの作成
    ax2.plot(xs_rel[:, 1], xs_rel[:, 2])
    ax2.set_xlabel("y")
    ax2.set_ylabel("z")
    ax2.set_xlim([-0.2, 0.2])
    ax2.set_ylim([-0.2, 0.2])

    ax3 = fig.add_subplot(133, aspect=1)  # 軌道のzx平面投影グラフの作成
    ax3.plot(xs_rel[:, 0], xs_rel[:, 2])
    ax3.set_xlabel("x")
    ax3.set_ylabel("z")
    ax3.set_xlim([-0.2, 0.2])
    ax3.set_ylim([-0.2, 0.2])

    plt.show()
    """


if __name__ == '__main__':
    cal_and_plot(1)
