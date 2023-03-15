# Determination-of-running-events
This program allows to determine the flight and contact times of a runner from a .csv file of the antero-posterior acceleration of the torso as well as its angular velocity in relation to this same axis.
This programme is based on the data received by an inertial unit placed on the runner's chest. The tests were carried out with a sensor from Movesens.

Please note: 
- this programme has been tested on people who run with a heel strike. It does not work on people running on their toes.
- The calculations do not start at the beginning of the recording. In fact, this delay allows the calculation to begin when the stride is correctly paced and the signal is clear (allowing a more precise calculation). If you want to change the timing of the start of the calculation, you must do so via the code. No interface has been created for this.
- The lower the running speed, the more inaccurate the calculation (because the acceleration curves are less pronounced). The program has been validated for speeds above 10km/h

Guide for use

1. You run with your inertial measurement system placed on your chest and make a recording.
2. You retrieve the recording files of the antero-posterior acceleration (axis crossing the torso from the back to the front) as well as the file of the angular acceleration around this same axis.
3. You import these files into your Matlab library
4. You launch the program, a dialog box will ask you to enter the names of the files corresponding to the acceleration and the angular velocity that you wish to process.
5. A graphical interface will appear. You will be able to access a button allowing you to look at the acceleration graph on which the characteristic points are placed allowing the calculation of the flight and contact time of the runner. Another button will bring up a table giving you these values.
