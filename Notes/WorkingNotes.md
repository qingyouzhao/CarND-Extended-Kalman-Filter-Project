

# Extended Kalman Filter 

![1555710226042](.\Kalman Filter design)
$$
\text{prediction} \\
x' = Fx+u \\
P'=FPF^T \\\
\\
\text{measurement update}\\
y=z-Hx \\
S=HPH^T+R \\
K=PH^TS^{-1} \\
x'=x+Ky \\
P'=(I-KH)P
$$

### Extension

$$
x' = f(x,u) \ \ \ where \ \ \ u =0 \\
y = z-h(x') using \ H_j instead\ of \ H \\
y = z_{radar} - h(x')
$$








## Terminologies

**Process noise** : $Q$ the uncertainty in change of velocity
$$
P' = FPF^T + Q \\
Q = GQ_vG^T\\
Q_v=
\begin{pmatrix} 
\sigma^2_{ax} & \sigma_{axy} \\
\sigma_{axy} & \sigma^2_{ay} \\
\end{pmatrix}
=
\begin{pmatrix} 
\sigma^2_{ax} 	& 0 			\\
0 				& \sigma^2_{ay} \\
\end{pmatrix}
\\
Q=
\begin{pmatrix} 
\frac{\Delta t^4} 4 \sigma^2_{ax} & 0 & \frac{\Delta t^3} 2 \sigma_{ax}^2 & 0 \\
0 & \frac{\Delta t^4} 4 \sigma^2_{ay} & 0 & \frac{\Delta t^3} 2 \sigma_{ay}^2\\
\frac{\Delta t^3} 2 \sigma^2_{ax} & 0 & {\Delta t^2} \sigma_{ax}^2 & 0 \\
0 & \frac{\Delta t^3} 2 \sigma^2_{ay} & 0 & {\Delta t^2}\sigma_{ay}^2 \\
\end{pmatrix}
=
\begin{pmatrix} 
\sigma^2_{ax} 	& 0 			\\
0 				& \sigma^2_{ay} \\
\end{pmatrix}
$$


Q

**Measurement noise**: uncertainty in sensor measurements

### Laser measurement

Laser measures a position $p_x$ and $p_y$

 $R$ represents the uncertainty in the position measurements we receive from the laser sensor.
$$
R=E[\omega\omega^T]=
\begin{pmatrix}
\sigma_{px}^2 & 0 \\
0			  & \sigma_{py}^2 \\
\end{pmatrix} \\
$$


### Radar measurement

$h$ function
$$
\text{h function specify how the predicted position and speed get mapped to the polar cordicnates of range, bearing and raage rate} \\
h(x')=
\begin{pmatrix}
\rho \\
\phi \\
\dot\rho \\
\end{pmatrix}
=
\begin{pmatrix}
\sqrt{{\rho'_x}^{2} + {\rho'_y}^{2}  } \\
\arctan(\frac{\rho'_y}{\rho'_x}) \\
\frac{\rho'_xv'_x+\rho'_yv'_y}{\sqrt{{\rho'_x}^{2} + {\rho'_y}^{2}  } } \\
\end{pmatrix}

\\
\text{}
y=z-h(x')
$$


$\rho$ is the distance to the pedestrian

$\phi$ is referenced counter clockwise x-axis is the bearing

$\dot\rho$ is the projection of the velocity $v $ onto the line $L$

We need a multi-dimensional Taylor series to approximate our $h$ function
$$
T(x) = f(a) + (x-a)^TDf(a) + \frac{1}{2!}(x-a)^TD^2f(a)*(x-a) + ...
$$
Similar to the 1-d Taylor series expansion
$$
T(x) = f(a) + \frac{f'(a)}{1!}(x-a) + \frac{f''(a)}{2!}(x-a)^2 +...
$$
$Df(a)$ is a Jacobian matrix, $D^2f(a)$ is a Hessian matrix

We only need the Jacobian because we are under the assumption $(x-a)$ is small.



Now we have derived Jacombian matrix $H_j$ for $h(x')$
$$
H_j = 
\begin{bmatrix}
\frac{p_x}{\sqrt{p_x^2+p_y^2}} &\frac {p_x}{\sqrt{p_x^2+p_y^2}} & 0 & 0\\
-\frac{p_y}{p_x^2+p_y^2} & \frac{p_x}{p_x^2+p_y^2}& 0 & 0\\
\frac{p_y(v_xp_y-v_yp_x)}{(p_x^2+p_y^2)^{3/2}} &\frac{p_x(v_yp_x-v_xp_y)}{(p_x^2+p_y^2)^{3/2}} &\frac{p_x}{\sqrt{p_x^2+p_y^2}} &\frac {p_x}{\sqrt{p_x^2+p_y^2}} \\
\end{bmatrix}
$$




## Code implementation snippets

```python
# Implement the filter function below
x = matrix([[0.], [0.]]) # initial state (location and velocity)
P = matrix([[1000., 0.], [0., 1000.]]) # initial uncertainty
u = matrix([[0.], [0.]]) # external motion
F = matrix([[1., 1.], [0, 1.]]) # next state function
H = matrix([[1., 0.]]) # measurement function
R = matrix([[1.]]) # measurement uncertainty
I = matrix([[1., 0.], [0., 1.]]) # identity matrix

// Kalman Filter variables, C++
VectorXd x;	// object state
MatrixXd P;	// object covariance matrix
VectorXd u;	// external motion
MatrixXd F; // state transition matrix
MatrixXd H;	// measurement matrix
MatrixXd R;	// measurement covariance matrix
MatrixXd I; // Identity matrix
MatrixXd Q;	// process covariance matrix



def kalman_filter(x, P):
    for n in range(len(measurements)):

        # measurement update
        z = matrix([[measurements[n]]])
        y = z - (H * x)
        S = H * P * H.transpose() + R
        K = P * H.transpose() * S.inverse()
        x = x + K*y
        P=(I-K*H)*P
        # prediction
        x = F*x + u
        P = F*P*F.transpose()
        
    return x,P
```



