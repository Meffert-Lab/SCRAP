# Setting up your platform prior to installation

SCRAP is installed through Conda, and code is cloned through Git.
While many users will already have Conda and Git, these additional instructions are provided for convenience.

## Windows
You will first install Windows Subsystem for Linux using PowerShell.
Run PowerShell as an administrator and then run the following command:

```
wsl --install
```

WSL will then install. After restart, a console will open. You may need to press "Enter" to resume the installation.
After installation completes, you will have to set a username and password.
After launching WSL, you can then install Git and Miniconda as follows.

```
sudo apt install git
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda init bash
```

You may need to close and reopen WSL following Miniconda installation.

## MacOS
### Using Homebrew
Homebrew is a great package manager from MacOS. You can install brew as follows:

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

*Note: read your terminal output carefully after Homebrew installation is complete. You may need to run additional commands to finish installation.*

Once brew is installed, you can install Git and Miniconda as follows.

```
brew install git
brew install Miniconda
```

You will have to initialize Conda for the shell you use. In most cases this will primarily be `zsh`, or alternatively `bash`.
You can check which shell you are using by launching a Terminal window and reading the title bar.

```
conda init zsh	#For zsh
conda init bash	#For bash
```

You may need to close and reopen your Terminal following Miniconda installation.

### Not using Homebrew
If you prefer not to use Homebrew, you can install Git as follows.

```
xcode-select --install
```

Miniconda installation differs for M1 Macs and non-M1 Macs.
#### Non-M1

```
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh --output Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

You will have to initialize Conda for the shell you use. In most cases this will primarily be `zsh`, or alternatively `bash`.
You can check which shell you are using by launching a Terminal window and reading the title bar.

```
conda init zsh  #For zsh
conda init bash #For bash
```

You may need to close and reopen your Terminal following Miniconda installation.

#### M1

```
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh --output Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

You will have to initialize Conda for the shell you use. In most cases this will primarily be `zsh`, or alternatively `bash`.
You can check which shell you are using by launching a Terminal window and reading the title bar.

```
conda init zsh  #For zsh
conda init bash #For bash
```

You may need to close and reopen your Terminal following Miniconda installation.

## Ubuntu
You can install Git and Miniconda as follows.

```
sudo apt install git
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

You may need to close and reopen your shell following Miniconda installation.
