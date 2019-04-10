import numpy as np

BJDo = 2458581.5 #08.04.2018
BJDf = 2458649.5 #15.06.2019

# MODIFY ONLY THIS PART!! ############################################ 
# Enter your To's, period, and output file name. This should be the
# name of your star without white spaces and should end with .BJD.
######################################################################

T0_file = open('t0s.txt','r')
T0s = T0_file.readlines()
per_file = open('periods.txt','r')
periods = per_file.readlines()
bjd_file = open('file_input_BJD','r')
file_names = bjd_file.readlines()

for i in range(len(T0s)):
    To = float(T0s[i].strip('\n'))
    per = float(periods[i].strip('\n'))
    file_out = file_names[i].strip('\n')

######################################################################
# Don't change the rest.
    To_up = To + round( (BJDo - To)/per)*per - per

    f = open(file_out, 'w')

    for j in range(1000):
        To_int = To_up + per*j
        if (To_int > BJDo):
            if (To_int < BJDf):
                f.write(str(j)+'   '+str(To_int)+'   \n')


    f.close()

