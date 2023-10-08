Sep 26
1. cntl for BE calculation;
2. control run for new domain;
3. radar obs producing;
4. radar da;
5. plot_reflectivity, etc.

Sep 28
1. make the ob.radar file, coding is done, bacth operation done, needed to be verified;
2. get the quality control instructions, need to test;

Hodiday Lists
1. shopping: glasses; dresses; 
2. the national postdoc fund form need to be filled QAQ;
3. literature review: (1)EnVar paper; (2) the AI for backgroud error convariance; (3) two science and nature paper;
4. get the make_obs package more readable and sent to Tao for code review;
5. Tornado fund done;
6. test the radar DA control run;
7. do some diagnose for the experiments.
8. Have fun with family and enjoy the rest.
9. Finish reading 钢铁是怎样炼成的。

Sep 30
1. make the ob.radar file, with the non-precipitation checking and velocity dealias (half done);
2. BG code review as below.

## note for the BG code review:
1. Line 308: how does kds define?是并行网格中的起始点嘛？
2. Line 317: should that be max(z_index, 1) instead of 0? should write comment for the hydro_mean.dat, how does the refl. threshold define.
3. Line 267: why 40? if the user's domain is 100, our code still work? does that mean the radar data is not available above level 40? 
4. Line 333: only take graupel into consideration, if they use the microphysical scheme that contains hail?
5. L359: count_rr > 1, maybe > 10 at least? or 100?
6. L327: 这个ze_qr的命名好奇怪哦，一眼看上去像是表示qr的等效反射率，可是实际上表示的是qr的数值的和，要不叫qr_sum?
7. L355: rootproc是什么意思？是输出多一点诊断信息的意思？
8. L546: it 是什么意思？
9. L552 不理解；
10. L476: radar_non_precip_rf 是多少呀 在哪里定义的？-888?-999?一般设置多少是非气象回波？有时只是缺测呢而不是非气象回波呢。
11. h按照500米一层我感觉可能有问题。
12. L603函数名不能理解。
13. L620就这个大于5度graupel不存在，那对于hail如果降到地面怎么办。
14. L682是我们加的嘛？那注册表里加了嘛？要说明一下单位呢？

Oct 6
1. BG的方案已提交；

Oct. 7
1. make the obs. especially finish the QC. document the process into slides. 
- (a)"best dealiasing method' and 'missing fourdd function':
    https://groups.google.com/g/pyart-users/c/o8oAHtEjdbE
- (b) Dealias doppler velocities using the Region Based Algorithm (easiest)
    https://arm-doe.github.io/pyart/examples/correct/plot_dealias.html#sphx-glr-examples-correct-plot-dealias-py
- (c) Dealiased, Cleaned Velocity Field (add the velocity texture, more reasonable)
    https://github.com/openradar/ams-open-radar-2023/blob/main/notebooks/pyart/pyart-corrections.ipynb

Oct. 8
1. radar obs in wrf format done, hand into Tao for Dart DA.
2. analyze the experiments. plot the refl. for the whole typhoon process.
3. try to grid the radar obs. into wrf grid, which would be better for visualization for comparison. 
