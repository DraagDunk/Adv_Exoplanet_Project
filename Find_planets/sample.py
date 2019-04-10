import numpy as np

star_nr = 1
RA_max = 24. ; RA_min = 0.
DEC_max = -85 ; DEC_min = -35
To_min = 2455000. ; To_max = 2457000.
Per_min = 5. ; Per_max = 15.
Tdur_min = 4.5 ; Tdur_max = 8.
mag_min = 8. ; mag_max = 10.

To_start = 2458178.5 # March 1st, 2018
To_end = 2458543.5 # March 1st, 2019

for k in range(star_nr):
    #star_name_file = 'StarNr_'+str(k+1)+'.dat'
    #star_name = 'StarNr_'+str(k+1)
    #BJDs = 'StarNr_'+str(k+1)+'.BJD'
    #RA = np.random.uniform(RA_min, RA_max)
    #DEC = np.random.uniform(DEC_max, DEC_min)
    #mag = np.random.uniform(mag_min, mag_max)
    #To = np.random.uniform(To_min, To_max)
    #Per = np.random.uniform(Per_min, Per_max)
    #Tdur = np.random.uniform(Tdur_min, Tdur_max)
    BJDs = 'AV_Dor.BJD'
    To = 
    Per = 
    epoch_nr = int((To_end - To_start)/Per)
    epoch_0 = int((To_start - To)/Per) + 1
    
    #f = open(star_name_file,'w')
    #f.write("%s %s %f %f %f %f %s\n" % ('1', star_name, Tdur, RA, DEC, mag, BJDs))
    #f.close()

    epochs = np.zeros(epoch_nr)
    Tos = np.zeros(epoch_nr)
    for j in range(epoch_nr):
        Tos[j] = To + epoch_0*Per + Per*j
        epochs[j] = j
    np.savetxt(BJDs, zip(epochs, Tos), fmt='%f')
    

    #print star_name, RA, DEC, Per, epoch_nr
