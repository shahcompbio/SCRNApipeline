from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.
additional_packages = ["appdirs","packaging","numpy","pypeliner","glob","os","argparse","rpy2"]
buildOptions = dict(packages = [], excludes = [])

base = 'Console'

executables = [
    Executable('pipeline.py', base=base)
]
additional_packages = []
setup(name='scrna_pipeline',
      version = '1.0',
      description = 'Single Cell RNA-Seq Pipeline',
      options = dict(build_exe = buildOptions),
      executables = executables)
