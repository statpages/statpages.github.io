# StatPages
Version control for John C Pezzullo's StatPages.info site


WORKING with GitHub

Each collaborator must have a GitHub account and a local running Git version (Windows, Linux, Mac)

1) WINDOWS INSTALL:

    http://windows.github.com

includes both command line tool and GUI
or

    https://git-scm.com/download/win
    
this installed effortlessly on my old XP. After setting the PATH, git run nicely from the command line


After install run the following commands:

1)

    git config --global user.name "John Doe"

2)

    git config --global user.email johndoe@example.com


That set You are ready for cloning. I recon that You are running a local internet server for testing before deploying (i.e. You do not edit directly on the server)?
If, cd to Your document root, else go (create) eg. C:/Documents/www

then run

1)

    git clone https://github.com/statpages/statpages.github.io.git

This will download the repostory to your PC. It will take some time because of the size of the repostory (Bill has a lot of stuff) and You are ready for work :-)

2)

    cd statpages.github.io

3)

    git status 

tells You that everything fine :-)



Commond commands:

1)
    git status

tells if anything changed

2)
    git add .

add all changes for commit

3)

    git commit -a -m 'comment on the particular commit'

commits/confirms Your edits for Git

4)

    git push

sends changes to GitHub. You will be prompted for Your Git password


It's essential that, BEFORE starting editing, You run a git pull to be sure You have the latest changes.
Although collaborators rarely will work on the exact same problem, conflicts may occur if You forget this step

Please contact me (soren@merser@gmail.com) in case You are in doubt, BEFORE commit/push