from azure.common.client_factory import get_client_from_cli_profile
from azure.mgmt.compute import ComputeManagementClient
from azure.mgmt.network import NetworkManagementClient
import os, uuid, sys
from azure.storage.blob import BlockBlobService, PublicAccess
from utils.config import *
import subprocess
import json
import tarfile
import shutil
import gzip


os.environ["AZURE_STORAGE_ACCOUNT"] = "singlecelldata"
os.environ["AZURE_STORAGE_KEY"] = "436b89a7-3b73-4644-a97b-949c4d0f19f5"

class ResultsInterface(object):

    def __init__(self, version="", dir="./tmp"):
        self.dir = dir + version
        try:
            os.makedirs(dir)
        except Exception as e:
            pass
        self.results = []
        self.version = version
        cmd = "az storage blob list --container-name cellranger{} --account-name scrnadata --auth-mode login".format(version)
        json_dump = subprocess.check_output(cmd.split())
        blobs = json.loads(json_dump)
        for blob in blobs:
            result = DataStorage(blob["name"])
            result.download_tenx(path=self.dir, version=version)
            self.results.append(result)

class TenxDataStorage(object):

    def __init__(self, sampleid, version="v3"):
        self.sampleid = sampleid
        self.storage_account = "scrnadata"
        self.container = "cellranger{}".format(version)
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token='?sv=2018-03-28&ss=bfqt&srt=sco&sp=rwdlacup&se=2021-03-19T02:52:48Z&st=2019-02-22T19:52:48Z&spr=https&sig=4oAGvIyqi9LPY89t21iWWp4XbbIgpmaldgcwWUOuf14%3D')
        self.tenx_path = None
        self.cache = ".cache"
        try:
            os.makedirs(self.cache)
        except Exception as e:
            pass

    def upload_cellranger(self, tenx):
        bam_tarball = tenx.bam_tarball()
        bam_tarball_name = os.path.split(bam_tarball)[-1]
        outs_tarball = tenx.outs_tarball()
        outs_tarball_name = os.path.split(outs_tarball)[-1]
        bam_blob_path = "{sampleid}/{local}".format(self.sampleid, bam_tarball_name)
        self.upload("bams" ,bam_blob_path, bam_tarball)
        outs_blob_path = "{sampleid}/{local}".format(self.sampleid, bam_tarball_name)
        self.upload("cellranger", outs_blob_path, bam_tarball)


    def unpack(self,path):
        tar = tarfile.open(path)
        tar.extractall()
        tar.close()
        os.rename("./outs",os.path.join(self.cache,self.sampleid))



    def download(self):
        local = ".cache/{}.tar.gz".format(self.sampleid)
        gzipped = "{}.tar.gz".format(self.sampleid)
        if not os.path.exists(local):
            self.block_blob_service.get_blob_to_path(self.container, gzipped, local)
            self.unpack(local)
        self.tenx_path = os.path.join(self.cache,self.sampleid)
        return self.tenx_path

    def upload(self,container, blob, local):
        print ("Uploading {} to {} in {}".format(container,local,blob))
        cmd = "az storage blob upload --container-name {} --file {} --name {} --auth-mode login".format(container,local,blob)
        subprocess.call(cmd.split())

class FastqDataStorage(object):

    def __init__(self, sampleid):
        self.sampleid = sampleid
        self.storage_account = "scrnadata"
        self.container = "rnaseq"
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token='?sv=2018-03-28&ss=bfqt&srt=sco&sp=rwdlacup&se=2021-03-19T02:52:48Z&st=2019-02-22T19:52:48Z&spr=https&sig=4oAGvIyqi9LPY89t21iWWp4XbbIgpmaldgcwWUOuf14%3D')
        self.tenx_path = None
        self.cache = ".cache"
        try:
            os.makedirs(self.cache)
        except Exception as e:
            pass

    def download_fastqs(self):
        fastqs = self.block_blob_service.list_blobs(self.container)
        for fastq in fastqs:
            print (fastq.name)
        # cmd = "az storage blob list --container-name  --account-name scrnadata --auth-mode login".format(self.fastq_container)
        # json_dump = subprocess.check_output(cmd.split())
        # blobs = json.loads(json_dump)
        # for blob in blobs:
        #     if blob["name"].startswith(self.sampleid):
        #         print("Downloading {}...".format(blob["name"]))
        #         blob_name = os.path.split(blob["name"])[1]
        #         local = os.path.join(self.datapath,blob_name)
        #         if not os.path.exists(local):
        #             cmd = "az storage blob download -n {blob} -f {path} --account-name scrnadata --container-name fastqs --auth-mode login".format(blob=blob_name, path=local)
        #             subprocess.call(cmd.split())
        # return self.datapath



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
    test_fastqs = FastqDataStorage("test-fastqs")
    test_fastqs.download_fastqs()
