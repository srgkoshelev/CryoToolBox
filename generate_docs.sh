echo 'Generating docs'
python3 -m pdoc --html CryoToolBox -o docs -c latex_math=True
mv docs/CryoToolBox/*.html docs
rm -rv docs/CryoToolBox
echo 'Docs generated'
while true; do
    read -p 'Commit and upload the docs? [y/n]' yn
    case $yn in
        [Yy]* )
            git add -A
            git commit -m 'Generated docs.'
            git push
            echo 'Committed and pushed upstream.'
            exit;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
