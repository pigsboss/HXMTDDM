## 介绍
HXMTDDM是硬X射线调制望远镜的配套成像软件包。
该软件包实现了加速直接解调方法，包含核心解调过程以及一系列功能性过程。

## 功能
目前支持的功能有：巡天演示程序、巡天数据处理、巡天图像查看、巡天图像保存。

### 巡天演示程序
为了便于我们进行软件包关键过程与模块的测试，开发了此演示程序。
巡天演示程序调用卫星轨道模拟器，生成巡天观测中的卫星姿态数据；
调用观测数据模拟器，生成模拟观测数据；
调用解调过程，对模拟观测数据进行解调，通过少量迭代重建较为粗糙的图像。

### 巡天数据处理
有MATLAB环境中交互式、有图形输出的巡天数据处理代码，以及编译过的可执行程序两种方式。

交互式代码的优点是便于调试，支持输入参数的可视化选取。

编译后的程序没有图形界面输出，不支持图形用户界面，仅有命令行环境的文字输出和文件读写。
适合使用预先定义的参数集对大量观测数据进行长时间批处理。

### 巡天图像查看与保存
在MATLAB桌面环境中，可以很方便的查看和保存重建过程中的图像以及其他类型的数据。
为了在其他环境，如Shell或通过网页脚本建立的会话等，方便的查看以及输出重要结果，我们开发了这个功能。
使用该功能，可以将用户指定天区指定迭代步骤的重建结果以图像及等高线图的形式显示在图形界面中，或者输出为EPS图像文件。


***

## Introduction
HXMTDDM is the imaging software package for Hard X-ray Modulation Telescope (HXMT), the first Chinese space telescope mission. This package implements an accelerated direct demodulation method (DDM) with the core demodulation procedure and a series of utility procedures.

## Functions
### All-sky survey demo
### All-sky survey worker
### All-sky survey viewer
### All-sky survey image printing
