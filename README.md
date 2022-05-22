# Solution pockets

## Run interactively with binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/georgkaufmann/Thesis_Gruber_SolutionPockets.git/main?filepath=index.ipynb)

## Versioning

We use [Git](https://git-scm.com/) for versioning.

Merge remote branch to local branch with the following steps:

1. Check local and remote branches:

~~~bash
$ git branch
* main
$ git branch -r
  github/main
~~~

2. Fetch status information for local and remote branches

~~~bash
$ git fetch github main
Username for 'https://github.com': XXXX
Password for 'https://georg.kaufmann@fu-berlin.de@github.com': 
From https://github.com/georgkaufmann/CO2
 * branch            main       -> FETCH_HEAD
~~~

3. List differences:

~~~bash
$ git log --oneline main..github/main
~~~

4. Merge remote branch to local branch:

~~~bash
$ git merge github/main
~~~

## Authors

* **Georg Kaufmann** - *Initial work* - [Georg Kaufmann](http://userpage.fu-berlin.de/~geodyn)

![](images/fu-logo.jpg)


## License

This project is licensed for classroom use only.

## Acknowledgments
