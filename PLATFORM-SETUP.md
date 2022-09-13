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
```

## MacOS
### Using Homebrew
Homebrew is a great package manager from MacOS. You can install brew as follows:

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Once brew is installed, you can install Git and Miniconda as follows.

```
brew install git
brew install Miniconda
```

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

#### M1

```
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh --output Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

## Ubuntu
You can install Git and Miniconda as follows.

```
sudo apt install git
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
