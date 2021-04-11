# Pull a pull request

pull PR #22 and to a new local branch called doc_tony
```
git fetch origin pull/22/head:doc_tony
git checkout doc_tony
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

# Reset to 20 mins ago
```
git reset --hard master@{"20 minutes ago"}
```

https://stackoverflow.com/questions/1223354/undo-git-pull-how-to-bring-repos-to-old-state


# filter large files

https://www.deployhq.com/git/faqs/removing-large-files-from-git-history

```
git filter-branch --force --index-filter \
  'git rm --cached --ignore-unmatch path/to/ceo.jpg' \
  --prune-empty --tag-name-filter cat -- --all
```
