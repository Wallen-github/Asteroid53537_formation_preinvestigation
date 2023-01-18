# Asteroid53537_formation_preinvestigation
We want to verify the formation of 53537-503955 pair. 


## Version 5

Primary rotation period currently: $P_{r1} = 72.74~hour$;

Primary rotation period at seperation: $P_{r2} = jP_a$; 

Evolution Time from seperation about 307~1467 kyrs;

Ratio of primary and orbit rotation: $j = \omega_A/n = 2,3$;

Mutual Orbit Period $(P_a)$ at seperate time is a function of mutual distance $(R)$: $P_a = \sqrt{\frac{4\pi^2R^3}{Gm_A}}$;

Yorp effect adoptes the expression in [(Rossi et al 2009)](https://www.sciencedirect.com/science/article/pii/S0019103509001109): $\tau_Y = \dot{\omega}=\frac{B_s F_s r_A}{a_s^2 \sqrt{1-e_s^2} ma} C_{\mathrm{Y}}$, in which $B_s = 2/3,F_s = 1\times10^{17} ~ kg~m/s^2$, and $a_s = 2.45 ~ AU,e_s = 0.079, C_{\mathrm{Y}} = 0.025$;

Based on the Yorp torque and primary;s rotation period at current time and speration time, we can estimate the evolution time after seperation:
$$T(j,R) = \frac{\frac{2 \pi}{P_{r 2}}-\frac{2 \pi}{P_{r 1}}}{\tau_Y} = \frac{1}{j\tau_Y}\sqrt{\frac{Gm_A}{R^3}} - \frac{2\pi}{P_{r1}\tau_Y}$$

Considering the Yorp torque is a constant, the envolution time is a function of $j$ and $R$.

![SepTime](SepTime.png)

