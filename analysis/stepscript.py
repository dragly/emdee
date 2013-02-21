# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 13:40:52 2013

@author: svenni
"""

from pylibconfig import Config
import subprocess
import atexit

processes = []

configDir = "configs/heating-freezing/"

#configFiles = [configDir + "step1.cfg", configDir + "step2.cfg", configDir + "step3.cfg", configDir + "step4.cfg", "configs/step5.cfg", "configs/step6.cfg"]
configFiles = [configDir + "step1.cfg", configDir + "step2.cfg", configDir + "step3.cfg", configDir + "step4.cfg", configDir + "step5.cfg", configDir + "step6.cfg", configDir + "step7.cfg", configDir + "step8.cfg", configDir + "step9.cfg"]

application = "../src-stan-build-Desktop_Qt_5_0_1_GCC_64bit-Debug/molecular-dynamics"

for configFile in configFiles:
    p = subprocess.Popen([application, configFile])
    processes.append(p)
    p.wait()