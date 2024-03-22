# LogFile.py

import time

# Logging messages throughout CrypSplice #



# Write Arguments 
## creating log file for logged messages
def write_arguments(outDir, args):
	logfile = open(outDir+'log.txt','w')
	logfile.write('# *********** Arguments *********** \n')
	args_dict = vars(args)
	for k in args_dict.keys():
	    if k in ["c1","c2","samples"]:
	        logfile.write("\t"+k+' : ' + "\t".join(args_dict[k])+'\n')
	    elif k in ['h','help']:
	        pass
	    else:
	        logfile.write("\t"+k+' : ' + str(args_dict[k])+'\n')
	logfile.write('# *********** ********* ************** \n')
	logfile.close()



# Log Message
## logging progress messages after each command 
def log_message(logfile_path, message):
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    with open(logfile_path, "a") as fileObj:
        fileObj.write("# " + message + ":" + localdate + '   ' + localtime + ' \n')
        print("# " + message + " : " + localdate + '   ' + localtime + '\n')
