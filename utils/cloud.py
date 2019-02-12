from azure.common.client_factory import get_client_from_cli_profile
from azure.mgmt.compute import ComputeManagementClient
from azure.mgmt.network import NetworkManagementClient
import os, uuid, sys
from azure.storage.blob import BlockBlobService, PublicAccess
from config import *
import subprocess
import json

from interface.tenxanalysis import TenxAnalysis

os.environ["AZURE_STORAGE_ACCOUNT"] = "singlecelldata"
os.environ["AZURE_STORAGE_KEY"] = "436b89a7-3b73-4644-a97b-949c4d0f19f5"

class DataStorage(object):

    def __init__(self, sampleid):
        self.sampleid = sampleid
        self.storage_account = "scrnadata"
        self.fastq_container = "rnaseq"
        self.bcl_container = "bcl"
        self.results_container = "cellranger"
        self.bams_container = "bams"
        self.block_blob_service = BlockBlobService(account_name='scrnadata',
                                                   account_key='436b89a7-3b73-4644-a97b-949c4d0f19f5')
        login = "az login -u nickceglia@gmail.com -p Aphahb_666"
        subprocess.call(login.split())
        self.datapath = "/home/nceglia/data/{sampleid}/".format(sampleid=self.sampleid)
        self.jobpath = "/home/nceglia/jobs/{sampleid}/".format(sampleid=self.sampleid)
        try:
            os.makedirs(self.datapath)
        except Exception as e:
            pass
        try:
            os.makedirs(self.jobpath)
        except Exception as e:
            pass
        self.jobpath = os.path.join(self.jobpath, self.sampleid)
        try:
            os.makedirs(self.jobpath)
        except Exception as e:
            pass

    def upload_results(self, tenx):
        bam_tarball = tenx.bam_tarball()
        bam_tarball_name = os.path.split(bam_tarball)[-1]
        outs_tarball = tenx.outs_tarball()
        outs_tarball_name = os.path.split(outs_tarball)[-1]
        bam_blob_path = "{sampleid}/{local}".format(self.sampleid, bam_tarball_name)
        self.upload("bams" ,bam_blob_path, bam_tarball)
        outs_blob_path = "{sampleid}/{local}".format(self.sampleid, bam_tarball_name)
        self.upload("cellranger", outs_blob_path, bam_tarball)


    def download_fastqs(self):
        cmd = "az storage blob list --container-name {} --account-name scrnadata --auth-mode login".format(self.fastq_container)
        json_dump = subprocess.check_output(cmd.split())
        blobs = json.loads(json_dump)
        for blob in blobs:
            if blob["name"].startswith(self.sampleid):
                print("Downloading {}...".format(blob["name"]))
                blob_name = os.path.split(blob["name"])[1]
                local = os.path.join(self.datapath,blob_name)
                if not os.path.exists(local):
                    cmd = "az storage blob download -n {blob} -f {path} --account-name scrna --container-name fastqs --auth-mode login".format(blob=blob_name, path=local)
                    subprocess.call(cmd.split())
        return self.datapath


    def download_tenx(self):
        tenx = "{}.tar.gz".format(self.sampleid)
        local = os.path.join(self.jobpath, tenx)
        cmd = "az storage blob download -n {blob} -f {path} --account-name scrna --container-name cellranger --auth-mode login".format(blob=tenx, path=local)
        return local

    def upload(self,container, blob, local):
        print ("Uploading {} to {} in {}".format(container,local,blob))
        cmd = "az storage blob upload --container-name {} --file {} --name {} --auth-mode login".format(container,local,blob)
        subprocess.call(cmd.split())


class VirtualMachine(object):

    def __init__(self):
        login = "az login -u nickceglia@gmail.com -p Aphahb_666"
        subprocess.call(login.split())
        compute_client = get_client_from_cli_profile(ComputeManagementClient)
        network_client = get_client_from_cli_profile(NetworkManagementClient)
        nic = network_client.network_interfaces.get(GROUP_NAME, NIC_NAME)
        print('\nCreating Linux Virtual Machine')
        vm_parameters = create_vm_parameters(nic.id, VM_REFERENCE['linux'])
        print(vm_parameters)
        async_vm_creation = compute_client.virtual_machines.create_or_update(GROUP_NAME, VM_NAME, vm_parameters)
        async_vm_creation.wait()

    def run_job(self):
        pass

    def check_status(self):
        pass

if __name__ == '__main__':
    testdata = DataStorage("SC_723")
    tenx = TenxAnalysis("/home/nceglia/jobs/SC_723/SC_723/outs")
    tenx.finalize()
    testdata.upload_results(tenx)
