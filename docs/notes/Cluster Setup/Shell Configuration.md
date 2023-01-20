# Shell Configuration
The first thing we can do is to add some things to our `.bashrc` file, which is the configuration file for the shell in most linux systems.

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
