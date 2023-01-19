In this section we will learn to set up <a href="https://mamba.readthedocs.io/en/latest/index.htm" target="_blank">Mamba</a> which is a fast package manager that can be used to create virtual environments and install packages without root permisions.

# Installing Mamba
Go to <a href="https://github.com/conda-forge/miniforge#mambaforge" target="_blank">this site</a> and download the appropriate version of **Mambaforge** for you system (if you are setting it up on a cluster use the *Linux x86_64(amd64)* option). If you are installing on a cluster you can use the following command to download the install script.

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
```

One the installation script is downloaded to your system you can run the following in the same folder it was downloaded in to install mamba.

```bash
sh Mambaforge-Linux-x86_64.sh
```

Once this command is run it will ask you to scroll through a license agreement, and then it will ask you some basic questions about where you would like mamba to be installed. It will look something like this. I like to install it into my `.local` directory, but it can be installed anywhere you want.

```bash
Welcome to Mambaforge 22.9.0-
In order to continue the installation process, please review the license agreement.
Please, press ENTER to continue
>>>
...
...
Do you accept the license terms? [yes|no]
[no] >>> yes
...
...
Mambaforge will now be installed into this location: /home/dd/mambaforge
- Press ENTER to confirm the location
- Press CTRL-C to abort the installation
- Or specify a different location below
[/home/dd/mambaforge] >>> /home/dd/.local/mambaforge
```

After this just agree to all the next steps until this final question comes up, and answer yes to it.

```bash
Do you wish the installer to initialize Mambaforge
by running conda init? [yes|no]
[no] >>> yes
```

Once this is finished you can source your `.bashrc` to activate you mamba environment.

```bash
source ~/.bashrc
```

If this was succesfull you should have `(base)` in front of your command prompt.