# cpp_nav_filters
A collection of C-Plus-Plus INS filters to perform fusion of GNSS and INS nav solutions

## To-Do
- [ ] Loosely Coupled INS EKF
    - [ ] Full state mechanization needs to be tested
- [ ] Tightly Coupled INS EKF
- [ ] GPS Least Squares
    - [ ] Test on multipe files/multiple epochs
        - [ ] Works on entire static data set.
        - [ ] Test GpsLeastSquares on new and exciting datasets
    - [X] Return DOP
    - [X] Resilience Checks
        - [X] Test for sufficient satellites before attempting to calculate a solution
    - [ ] Weighting matrix
        - [ ] Use CN0 estimate
        - [ ] Use elevation angle
    - [ ] Set mask angle on satellites?
- [ ] Common 
    - [X] Struct of ephemeris for each satellite
        - [X] Struct contains vectors of variable dimension for each ephemeride
        - [X] A vector of same variable dimension that contains ID for each SV
    - [X] Resilience Checks
        - [X] Ensure common has ephemerides to calcualte SV PVT states for the sv_id satellite
    - [ ] General Navigation Utilities
        - [ ] Accelerometer Levelling (Tested on single data set and matches MATLAB results)
        - [ ] Somigliana Gravity Model (groves p. 72)
        - [ ] Rotation matrices to euler angles and vis. versa. 
- [ ] Plotting in CPP
    - [ ] Currently Using matplotplusplus as plotting tool in c++

## Working Inputs and Outputs of Various bits of Cpp Nav Filt
- Common
    - receiveSvEphem
        - Inputs: vector of ephemerides and int of satellite (0->31)
        - Outputs: None
        - Purpose: Common contains a matrix of ephemerides for all GPS Sv's, this function serves to store all ephemerides in one place 
    - sendSvEphem
        - Inputs: ephem_out and desired sv (0->31)
        - Outputs: ephem_out
        - Purpose: If a user requires SV ephemeris, this function will grab the correct ephemeris from the common struct and pass by reference it to the user 
    - sendUnitVectors
        - Inputs: state estimate (PVT), SVPVT, H
        - Outputs: H is the matrix of unit vectors (pass by reference)
        - Purpose: Compute a set of unit vectors to N satellites given a user PVT estimate and SV PVT states
    - sendMeasEst
        - Inputs:  state estimate (PVT), SVPVT, Y_hat
        - Outputs: Y_hat (pass by reference)
        - Purpose: Given a user PVT state estimate and SVPVT states, compute estiamted pseudoranges and pseudorange rates
    - sendSvStates
        - Inputs: sv desired (0->31), transmit time, transit time
        - Outputs: SV PVT (x,y,z,x_dot,y_dot,z_dot,clk correction)
        - Purpse: Compute the ECEF position and velocity and clock correction terms for a satellite at a specific point in time
      
- Gps Least Squares
    -  sendStateEstimate
        - Inputs: Meaurement Vector (pseudoranges and pseudorange rates), SVPVT matrix, common class, state estimate
        - Outputs: (x,y,z,b,x_dot,y_dot,z_dot,b_dot) user rcvr PVT states
        - Purpose: Use Newton-Raphson to compute a user PVT solution at a single epoch from a given set of pseudorange and pseudorange rate measurements and a matrix of SV pvt states 

## Test Data
Test data is pulled from Auburn University Fundamentals of GPS Navigation lab data as it is the data I am most familiar with:
[Fundamentals of GPS Navigation Class Website](https://eng.auburn.edu/~dmbevly/fund_gps/)
