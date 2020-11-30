'''
Description: Run CellRanger Count
Author: Liuzj
Date: 2020-10-12 11:08:25
LastEditTime: 2020-11-27 11:59:36
LastEditors: Liuzj
'''
import yaml
import sh
from loguru import logger

def parseCellRangerParameters(parameterPath):
    """
    need additional parameters:
        cellRangerPath
        outDir
    """
    with open(parameterPath) as fh:
        parameters = yaml.load(fh)
        shNeedParasDict = {}
        shNeedParasDict['cellRangerPath'] = parameters.pop('cellRangerPath')
        shNeedParasDict['cellRangerOutDir'] = parameters.pop('outDir')
        shParasList = []
        for crOption, crValue in parameters.items(): #cr cellRanger count
            if crValue == '':
                shParasList.append(f'--{crOption}')
            else:
                shParasList.append(f'--{crOption}={crValue}')
        shNeedParasDict['cellRangerParas'] = shParasList
    return shNeedParasDict

def runCellRangerParameters(shNeedParasDict):
    try:
        sh.mkdir(shNeedParasDict['cellRangerOutDir'])
    except:
        logger.warning('{shNeedParasDict["cellRangerOutDir"]} existed !')
    sh.cd(shNeedParasDict['cellRangerOutDir'])
    sh.Command(shNeedParasDict['cellRangerPath']).count(shNeedParasDict['cellRangerParas'])

def runCellRanger(parameterPath):
    """
    run cellRanger count
    need additional parameters:
        cellRangerPath
        outDir
    """
    shNeedParasDict = parseCellRangerParameters(parameterPath)
    logger.info('cellRanger run Start')
    runCellRangerParameters(shNeedParasDict)
