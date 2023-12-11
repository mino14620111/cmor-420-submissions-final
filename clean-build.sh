
mkdir -p .build

pdflatex -output-directory .build "$project.tex"

mv .build/"$project.pdf" 

git add clean-build.sh
git commit -m "Add clean-build script"
mkdir -p cmor-420-520-submissions/homework-1
mv clean-build.sh cmor-420-520-submissions/homework-1/



