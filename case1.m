%{ 
TOPIC:      Odometric Error Propagation.
NAME:       Naol Samuel Erega.
STUDENT ID: 20210889
%}

% Parameters
steps=200;
sigma_level=1;        % confidence level(1sigma= 0.6827)
ds= 0.2;              % robot mov(m)
dtheta=0;             % case I
b=0.2;                % base length(m)

% Initialize position and pose covariance
x = 0;
y = 0;
theta = 0;

Cp = [1^2 0 0;
      0 1^2 0;
      0 0 0];

% Covariance of wheel slip
Cw = [0.1^2, 0;
      0, 0.1^2];

% Setup single plot for trajectory and error ellipses
figure;
hold on; grid on;
title('dtheta 0 deg');
xlabel('X position');
ylabel('Y position');

% Loop through wheel movements
for I = 1:1:steps
    oldx = x;
    oldy = y;

    % Update position and heading
    [x, y, theta, Cp] = pos_update(x, y, theta, Cp,Cw,ds,dtheta, b);

    % Plot trajectory
    plot([oldx x], [oldy y], 'b--');
    if mod(I,10) == 0
        quiver(x, y, cos(theta), sin(theta), 'b');
    end
    
    if(rem(I,20)==0)
         % Calculate and plot error ellipses (for x-y plane only)
        R = chol(Cp(1:2, 1:2), 'lower');
        theta_vals = linspace(0, 2*pi, 100);
        unit_circle = [cos(theta_vals); sin(theta_vals)];
        ellipse_points = R * sigma_level * unit_circle;
        plot(x + ellipse_points(1,:), y + ellipse_points(2,:), 'r --');
    end
end
drawnow;
xticks(min(xlim):2:max(xlim));


function [x, y, theta, Cp] = pos_update(x, y, theta, Cp,Cw,ds,dtheta, b)
    % Calculate the average movement and rotation change
    dx = ds * cos(theta + dtheta / 2);
    dy = ds * sin(theta + dtheta / 2);

    % Partial derivatives of kinematics with respect to x, y, theta
    Pp = [1, 0, -ds * sin(theta + dtheta / 2);
          0, 1, ds * cos(theta + dtheta / 2);
          0, 0, 1];

    % Partial derivatives of kinematics with respect to wheel motion
    Pw = [0.5 * cos(theta + dtheta/2) - 0.5/b * ds * sin(theta + dtheta/2), ...
          0.5 * cos(theta + dtheta/2) + 0.5/b * ds * sin(theta + dtheta/2);
          0.5 * sin(theta + dtheta/2) + 0.5/b * ds * cos(theta + dtheta/2), ...
          0.5 * sin(theta + dtheta/2) - 0.5/b * ds * cos(theta + dtheta/2);
          1 / b, -1 / b];

    % Update position
    x = x + dx;
    y = y + dy;
    theta = theta + dtheta;

    % Update covariance matrix
    Cp = Pp * Cp * Pp' + Pw * Cw * Pw';
end
