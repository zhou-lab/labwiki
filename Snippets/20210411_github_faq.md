# Pull a pull request

pull PR #22 and to a new local branch called doc_tony
```
git pull origin pull/22/head:doc_tony
```

see
https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/checking-out-pull-requests-locally

# fetch and FETCH_HEAD
FETCH_HEAD is a temporary head
```
git pull origin master
```
is equal to
```
git fetch origin master
git merge FETCH_HEAD
```

# How to interpret refspec?

https://stackoverflow.com/questions/44333437/git-what-is-refspec

https://stackoverflow.com/questions/27567846/how-can-i-check-out-a-github-pull-request-with-git
