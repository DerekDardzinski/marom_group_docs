# Installing the Helix Text Editor
In most linux systems the default text editor is Vi/Vim/Nano. These are all either difficult to work in or difficult to configure. I would recommend installing the <a href="https://helix-editor.com/" target="_blank">Helix text editor</a> since it comes pre-installed with language-server functionalities and lots of nice themes.

I found that the most reliable way to install on our clusters is building it from source. The first step is to install `rust` using mamba

```bash
mamba install rust
```

Now you can make a new directory in your home folder called `pkgs`

```bash
mkdir ~/pkgs
```

and go into the directory via:

```bash
cd ~/pkgs
```

Then you can run the following command to get the source code for the helix text editor

```bash
git clone https://github.com/helix-editor/helix
```

```bash
cd helix
```

```bash
cargo install --locked --path helix-term
```

```bash
mkdir -p ~/.config/helix
```

```bash
ln -s $PWD/runtime ~/.config/helix/runtime
```

```bash
echo 'export PATH="${HOME}/.cargo/bin:$PATH"' >> ~/.bashrc
```

```bash
source ~/.bashrc
```

This should give you a working helix installation and the editor can be opened using the command `hx`

A reasonable starting configuration is the following, and it can be inserted in the config file via:

```bash
hx ~/.config/helix/config.toml
```

Copy and paste this into config.toml thn run `:wq` to save the config file.

```toml
theme = "catppuccin_mocha"

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

The last thing to do is install a python and bash language server so you can get autocomplete and linting and autoformatting:

```bash
mamba install python-lsp-server black nodejs
```

```bash
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
