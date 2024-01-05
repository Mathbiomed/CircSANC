--------------------------------------------

Matlab codes for a mathematical model for the circadian regulation of sinoatrial nodal cell pacemaking function

This model was developed using the Morotti et al. model as the baseline, as described in: 

Pan Li and Jae Kyoung Kim, Circadian regulation of sinoatrial nodal cell pacemaking function: dissecting the roles of autonomic control, body temperature, and local circadian rhythmicity, 2023

To use the model: 

SingleRun_ZT.m                  to simulate SANC pacemaking at a given Zeitgeber time (ZT)
PerformanceAmplitude_ZTloop.m   to calculate circadian (24h) amplitude and performance of SANC 
				pacemaking
Paraspace.m                     to calculate the Ca oscillator and membrane oscillator parameter 
				space map at a given ZT
mouse_SAM_eccODEfile.m          ODE file for adult mice
mouse_SAM_eccODEfile_aging.m    ODE file for aged mice   

--------------------------------------------

Matlab code of the Morotti et al. model of the sinoatrial node myocyte in mouse updated to investigate the role of the Dnajb6 gene in sick sinus syndrome

The Morotti et al. model (https://github.com/drgrandilab/Ding-et-al-2022-mouse-sinoatrial-model) was originally developed by modifying the Kharche et al. model (Am. J. Physiol. Heart Circ. Physiol. 2011, 301, H945–H963) to reproduce an experimental dataset collected in the Proenza Lab.


Please cite the following papers when using the Morotti et al. code:

Morotti S, Ni H, Peters CH, Rickert C, Asgari-Targhi A, Sato D, Glukhov AV, Proenza C, Grandi E. Intracellular Na+ Modulates Pacemaking Activity in Murine Sinoatrial Node Myocytes: An In Silico Analysis. Int J Mol Sci. 2021 May 26;22(11):5645. doi: 10.3390/ijms22115645.

Ding Y, Lang D, Yan J, Bu H, Li H, Jiao K, Yang J, Ni H, Morotti S, Le T, Clark KJ, Port J, Ekker SC, Cao H, Zhang Y, Wang J, Grandi E, Li Z, Shi Y, Li Y, Glukhov AV, Xu X. A phenotype-based forward genetic screen identifies Dnajb6 as a sick sinus syndrome gene. Elife. 2022 Oct 18;11:e77327. doi: 10.7554/eLife.77327.


License: MIT License.
