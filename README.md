# Asteroid53537_formation_preinvestigation
We want to verify the formation of 53537-503955 pair. 


## Version 5
现在主星自转周期：$P_{r1} = 72.74~hour$

分裂时间区间：307～1467 千年

分裂时主星/轨道频率比：$j = \omega_A/n = 2,3$

双星分裂时的轨道周期($P_a$)与分裂时刻的双星间距($R$)有关：$P_a = \sqrt{\frac{4\pi^2R^3}{Gm_A}}$

Yorp效应采用[(Rossi et al 2009)](https://www.sciencedirect.com/science/article/pii/S0019103509001109)中的表达式：$\tau_Y = \dot{\omega}=\frac{B_s F_s r_A}{a_s^2 \sqrt{1-e_s^2} ma} C_{\mathrm{Y}}$, 其中$B_s = 2/3,F_s = 1\times10^{17} ~ kg~m/s^2,a_s = 2.45~AU,e_s = 0.079,C_{\mathrm{Y}} = 0.025$

分裂时主星的自转周期：$P_{r2} = jP_a$

基于Yorp力矩和主星现在和分裂时的自转周期，我们可以推算出双星分裂后的演化时间：
$$
T(j,R) = \frac{\frac{2 \pi}{P_{r 2}}-\frac{2 \pi}{P_{r 1}}}{\tau_Y} = \frac{1}{j\tau_Y}\sqrt{\frac{Gm_A}{R^3}} - \frac{2\pi}{P_{r1}\tau_Y}
$$
在Yorp力矩为常数时，分裂后到此刻的演化时间是$(j,R)$的函数，即共振比和分裂时刻双星间距。
