

The first time you want to run the simulations after downloading the
code, execute

>>  gsimStartup

This command will create a file named gsim.m, which contains the
function that you should invoke to execute any experiment.

All experiments are located in the folder Experiments, and its 
name is "DecentralizedProjectionsExperimentsTSP.m". The file 
contains one function per experiment, which are named experiment_N,
where N is a number identifying the experiment. 
For example to execute experiment 3001, which plots Fig.1 of the paper,
type 

>>  gsim(0,3001)

You will see Fig.1 of the paper that has been displayed on the screen.