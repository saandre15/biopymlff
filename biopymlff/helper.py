from tensorflow.python.client import device_lib

def get_avaliable_gpu():    
    local_devices_protos = device_lib.list_local_devices()
    return [x.name for x in local_device_protos if x.device_type == 'GPU']