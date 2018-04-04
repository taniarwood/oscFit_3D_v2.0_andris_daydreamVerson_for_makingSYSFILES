#!/usr/bin/env python


def runScript(job_script, fit_name, job_range):
    
    # Introduce more ways in which this can be run


    from IPython.parallel import Client, interactive
    import iparallel
    rc = Client(profile = 'wgs3')
    lview = rc.load_balanced_view()

    def runMCFit(job_script, file_name, task_nr):
        import os
        file_name += str(int(task_nr))
        qsub_string = '  '.join(['python',job_script, file_name])
        print qsub_string
        os.system(qsub_string)


    result = lview.map_async(runMCFit, [job_script]*len(job_range), [fit_name]*len(job_range), job_range)
    iparallel.waitOn(result)


if __name__ == '__main__':

    fit_name   = 'MCTest_PPC_Vacuum_Cascades_StatFluct/MCtest_PPC_DirectHI_Cascades_NoGammaPrior_StatFluct'
    job_range  = range(23,101)
    job_script = 'fitMC_mod.py'

    # Make a copy of the job_script. Call it aux. Run it, and store a copy in the directory.

    runScript(job_script, fit_name, job_range)
