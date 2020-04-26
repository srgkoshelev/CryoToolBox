# git -C wiki pull
pdoc --html heat_transfer -o .
# rm -v heat_transfer/index.md
mkdir docs
mv -v heat_transfer/*.html docs
rm -rv heat_transfer
# for filename in wiki/*.md; do
#     pandoc -f markdown -t gfm "$filename" -o "$filename"
# done
git add -A
git commit
git -C docs push
