echo 'Generating docs'
pdoc --html heat_transfer -o docs
mv docs/heat_transfer/*.html docs
rm -rv docs/heat_transfer
echo 'Docs generated'
while true; do
    read -p 'Commit and upload the docs? [y/n]' yn
    case $yn in
        [Yy]* )
            git add -A
            git commit -m 'Generated docs.'
            git push
            echo 'Committed and pushed upstream.'
            ;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
