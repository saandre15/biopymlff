import multiprocessing
import subprocess

from tensorflow.python.client import device_lib


def get_avaliable_gpu():    
    local_devices_protos = device_lib.list_local_devices()
    return [x.name for x in local_device_protos if x.device_type == 'GPU']

def get_cpu_count():
    return multiprocessing.cpu_count()

def get_gpu_count():

    try:
        proc = subprocess.Popen("nvidia-smi --list-gpus | wc -l", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
    except OSError as err:
        # Actually this may never happen with shell=True, since
        # probably the shell launches successfully.  But we soon want
        # to allow calling the subprocess directly, and then this
        # distinction (failed to launch vs failed to run) is useful.
        msg = 'Failed to execute "{}"'.format(command)
        raise EnvironmentError(msg) from err
    # stderr, stdout = proc.communicate()
    errorcode = proc.wait()

    if errorcode:
        path = os.path.abspath(self.directory)
        msg = ('Calculator "{}" failed with command "{}" failed in '
                '{} with error code {}'.format(self.name, command,
                                                path, errorcode))
    return int(str(proc.stdout.readline(1)).replace("b", "").replace("'", ""))

def run_parallel(fn: list):
    """ 
    Takes a list of callables to be runned in parallel.

    Automated solution to the max job duration problem
    """
    pass

