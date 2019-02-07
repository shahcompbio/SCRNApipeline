from azure.common.client_factory import get_client_from_cli_profile
from azure.mgmt.compute import ComputeManagementClient
from azure.mgmt.network import NetworkManagementClient
from config import *

class SampleData(object):

    def __init__(self, sampleid):
        self.sampleid = sampleid
        self.storage_account = "scrnadata"
        self.fastq_container = "fastqs"
        self.bcl_container = "bcl"
        self.results_container = "cellranger"
        self.bams_container = "bams"
        self.client = get_client_from_cli_profile(ComputeManagementClient)

    def results(self):
        return None

class VirtualMachine(object):

    def __init__(self):
        compute_client = get_client_from_cli_profile(ComputeManagementClient)
        network_client = get_client_from_cli_profile(NetworkManagementClient)
        nic = network_client.network_interfaces.get(GROUP_NAME, NIC_NAME)
        print('\nCreating Linux Virtual Machine')
        vm_parameters = create_vm_parameters(nic.id, VM_REFERENCE['linux'])
        print(vm_parameters)
        async_vm_creation = compute_client.virtual_machines.create_or_update(
            GROUP_NAME, VM_NAME, vm_parameters)
        async_vm_creation.wait()

if __name__ == '__main__':
    testvm = VirtualMachine()
