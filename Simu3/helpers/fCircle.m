function H=fCircle(center,radius,varargin)
% FCIRCLE   Draws a circle in 2D or 3D.
%   H = FCIRCLE(CENTER,RADIUS) draws a circle in 2D. CENTER and RADIUS
%   specifies the properties of the circle. H is a handle to the object.
%
%   H = FCIRCLE(CENTER,RADIUS,NORMAL) draws a circle in 3D. NORMAL is the
%   normal vector for the circle (default is [0 0 1]').
%
%   H = FCIRCLE(CENTER,RADIUS,NORMAL,STYLE) specifies the style of the
%   object. Default is 'k-';
%
%   H = FCIRC(CENTER,RADIUS,NORMAL,STYLE,POINTS) specifies the number of
%   points to use in the circle. Default is 1000.

% Input argument check
if(nargin < 2 || nargin > 5),
    error('Wrong number of input arguments.');
elseif(nargin==2)
    n = [0 0 1]';
    style = 'k-';
    points = 1000;
elseif(nargin==3)
    n = varargin{1};
    style = varargin{2};    
    points = 1000;
else
    n = varargin{1};
    style = varargin{2};    
    points = varargin{3};
end

% Correct n
if(size(n,1) == 1 && size(n,2) == 3)
    n = n';
elseif(size(n,1) ~= 3 || size(n,2) ~= 1);
    error('Normal must be a 3x1 vector');
end


% Plotting points
theta=linspace(0,2*pi,points);

% Plot the circle
if(all(n == [0 0 1]')) % 2D-plot
    rho=ones(1,points)*radius;
    [X,Y] = pol2cart(theta,rho);
    X=X+center(1);
    Y=Y+center(2);
    H=plot(X,Y);
else % 3D-plot 
    % ---------------------------------------------------------------------
    % Find orth. vectors ; (n,v,w) is a ON-base ; n = normal from circle
    [thetaN, phiN, dummy] = cart2sph(n(1),n(2),n(3));
    
    % Rotates around z-axis, with the angle thetaN
    Qz=[cos(thetaN)     -sin(thetaN)    0 ;
        sin(thetaN)     cos(thetaN)     0 ;
        0                   0                   1];
    % Rotates around the y-axis with the angle phiN
    Qy=[cos(phiN)       0                   -sin(phiN) ;
        0              1                   0 ;
        sin(phiN)       0                   cos(phiN)];

    v = Qz*Qy*[0 1 0]';
    w = Qz*Qy*[0 0 1]';
    % ---------------------------------------------------------------------
    
    X = center(1) + radius*(v(1)*cos(theta) + w(1)*sin(theta));
    Y = center(2) + radius*(v(2)*cos(theta) + w(2)*sin(theta));
    Z = center(3) + radius*(v(3)*cos(theta) + w(3)*sin(theta));
    
    H = plot3(X,Y,Z,style);
end
    
    