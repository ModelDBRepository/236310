This model achieves electrotonic transform and computes mean inward and outward
attenuation from 0 to 500 Hz input and randomly activates synapses along
dendrites to simulate AMPAR mediated EPSCs.

For electrotonic analysis, in Elec folder, the entry file is
MSNelec_transform.hoc. 
For EPSC simulation, in Syn folder, the entry file is
randomepsc.hoc. 
Run read_EPSCsims_mdb_alone.m next with the simulated parameter
values specified to compute the mean EPSC.