# Cluster Setup
In this section we will learn how to set up a proper development environment on a linux cluster.

## Installing a Package Manager
In this section we will learn to set up <a href="https://mamba.readthedocs.io/en/latest/index.htm" target="_blank">Mamba</a> which is a fast package manager that can be used to create virtual environments and install packages without root permisions.

### Installing Mamba
Go to <a href="https://github.com/conda-forge/miniforge#mambaforge" target="_blank">this site</a> and download the appropriate version of **Mambaforge** for you system (if you are setting it up on a cluster use the *Linux x86_64(amd64)* option). If you are installing on a cluster you can use the following command to download the install script.

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
```

Once the installation script is downloaded to your system you can run the following in the same folder it was downloaded in to install mamba.

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

## Installing a Text Editor
### Helix Text Editor
In most linux systems the default text editor is Vi/Vim/Nano. These are all either difficult to work in or difficult to configure. I would recommend installing the <a href="https://helix-editor.com/" target="_blank">Helix text editor</a> since it comes pre-installed with language-server functionalities and lots of nice themes. If you aren't familiar with working in a vim-like editor try going through the tutorial after you install the editor `hx --tutor`.

#### Installing Rust
I found that the most reliable way to install on our clusters is building it from source. The first step is to install `rust` using mamba

```bash
mamba install rust
```

#### Downloading and Building the Source Code
Now you can make a new directory in your home folder called `pkgs` and `cd` into it.

```bash
mkdir ~/pkgs
cd ~/pkgs
```

Then you can run the following commands to get the source code for the helix text editor and instal it using cargo, the rust package manager.

```bash
git clone https://github.com/helix-editor/helix
cd helix
cargo install --locked --path helix-term
mkdir -p ~/.config/helix
ln -s $PWD/runtime ~/.config/helix/runtime
echo 'export PATH="${HOME}/.cargo/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

This should give you a working helix installation and the editor can be opened using the command `hx`

#### Configuring Helix
A reasonable starting configuration is the following, and it can be inserted in the config file via:

```bash
hx ~/.config/helix/config.toml
```

Copy and paste this into the config.toml file using `Ctrl-Shift-V`. Then you can run `:wq` to save and exit from the config file.

```toml
theme = "catppuccin_mocha" # there are many options, this is just the one I like to use. Try :theme <TAB> to see all available themes.

[editor]
true-color = true
line-number = "relative"
mouse = false

[editor.cursor-shape]
insert = "bar"
normal = "block"
select = "underline"

[editor.file-picker]
hidden = false

[keys.normal]
C-o = ":open ~/.config/helix/config.toml" # Maps the Ctrl-o to opening of the helix config file
g = { a = "code_action" } # Maps `ga` to show possible code actions
"ret" = ["open_below", "normal_mode"] # Maps the enter key to open_below then re-enter normal mode

[keys.insert]
"A-f" = "normal_mode" # Maps Alt-X to enter normal mode
```

Now you can reopen any file with `hx` and you should see a nice looking terminal-based text editor.

#### Adding Language-Servers and Auto-Formatters
The last thing to do is install a python and bash language server so you can get autocomplete and linting and autoformatting:

```bash
mamba install python-lsp-server black nodejs
npm i -g bash-language-server
```

For the python language settings you can add the following the the languages.toml file:

```bash
hx ~/.config/helix/languages.toml
```

```toml
[[language]]
name = "python"
auto-format = true
[language.formatter]
command = "black"
args = [ "-", "--quiet", "--line-length=79" ]
```

Now you should have a working text editor with language support for python and bash. 

**It is important to reinstall the `python-lsp-server` and `black` packages in every new mamba environment that you make in order for autocomplete to work for the python packages that are installed in that environment.**


## Configuring your Shell
### Adding comands to the .bashrc
The `.bashrc` file lives in your home directory and is the configuration file for you shell. Below are some commands that you can add to your `.bashrc` to make working in your shell a better experience than looking at plain white text in a terminal.

```bash
ulimit -s unlimited

if [ -f /etc/bashrc ]; then
		. /etc/bashrc
fi

# Add colors to the output of ls
export LS_OPTIONS='--color=auto'
eval "$(dircolors -b)"
alias ls='ls $LS_OPTIONS'

# Create a more informative propmt
PS1='\[\e[1;36m\]$(whoami)(\h) \[\e[1;35m\]${PWD#"${PWD%/*/*}/"}\[\e[00m\] ';

# Add your bin folder to your path
if [ ! -d "${HOME}/bin" ]; then
	mkdir "${HOME}/bin"
fi
export PATH="${HOME}/bin:$PATH"

if [ -f ~/.bash_aliases ]; then 
		. ~/.bash_aliases; 
fi
```

Once you add this, you can source the `.bashrc` file and you should see some color in your terminal.

```bash
source ~/.bashrc
```
