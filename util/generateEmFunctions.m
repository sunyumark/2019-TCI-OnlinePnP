function paramsOut = generateEmFunctions(p)
% Generate Green's functions

% meshgrid the pixel locations
[XPix, YPix] = meshgrid(p.x, p.y);

% handle to the Hankel function
hankFun = @(x) 1j*0.25*besselh(0, 1, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% locations of transmitters and receivers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% location of the transmissions
transmitterAngles = linspace(0, 359, p.numTrans)*pi/180; % [rad]
x_transmit = p.sensorRadius * cos(transmitterAngles); % [m]
y_transmit = p.sensorRadius * sin(transmitterAngles); % [m]

% location of the receivers
receiverAngles = linspace(0, 359, p.numRec)*pi/180; % [rad]
x_receive = p.sensorRadius * cos(receiverAngles);
y_receive = p.sensorRadius * sin(receiverAngles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance data between sensors and pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mask
p.receiverMask = true(p.numTrans, p.numRec);

% receiver to pixel
diff_x_rp = repmat(XPix, [1, 1, p.numRec])-...
    repmat(reshape(x_receive, [1, 1, p.numRec]), [p.Ny, p.Nx]);
diff_y_rp = repmat(YPix, [1, 1, p.numRec])-...
    repmat(reshape(y_receive, [1, 1, p.numRec]), [p.Ny, p.Nx]);
distanceRecToPix = sqrt(diff_x_rp.^2 + diff_y_rp.^2);

% transmitter to pixel
diff_x_tp = repmat(XPix, [1, 1, p.numTrans])-...
    repmat(reshape(x_transmit, [1, 1, p.numTrans]), [p.Ny, p.Nx]);
diff_y_tp = repmat(YPix, [1, 1, p.numTrans])-...
    repmat(reshape(y_transmit, [1, 1, p.numTrans]), [p.Ny, p.Nx]);
distanceTransToPix = sqrt(diff_x_tp.^2 + diff_y_tp.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.uincDom = hankFun(p.kb*distanceTransToPix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensor Green's functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sensorGreensFunction = hankFun(p.kb*distanceRecToPix);
p.sensorGreensFunction = (p.kb^2)*sensorGreensFunction;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain Green's functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Locations of the Green's function's pixels
xGreen = (-p.Nx:p.Nx-1)*p.dx;
yGreen = (-p.Ny:p.Ny-1)*p.dy;

% Meshgrid the pixel locations
[XGreen, YGreen] = meshgrid(xGreen, yGreen);
R = sqrt(XGreen.^2+YGreen.^2);

%%% generate Hankel function and remove singularity
domainGreensFunction = hankFun(p.kb*R);
domainGreensFunction(p.Ny+1,p.Nx+1) = quad2d(@(x,y) hankFun(p.kb*sqrt(x.^2+y.^2)),...
    -p.dx/2, p.dx/2, -p.dy/2, p.dy/2)/(p.dx*p.dy);
p.domainGreensFunction = (p.kb^2)*domainGreensFunction;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramsOut = p;