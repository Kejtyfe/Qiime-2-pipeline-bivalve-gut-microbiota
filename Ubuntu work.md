So first of all, we need to open our Bash - Ubuntu (Windows Subsystem for Linux)

Now activation of the environment  ദ്ദി◝ ⩊ ◜.ᐟ

```bash
conda activate qiime2-amplicon-2024.5 
```
Now we will create a file for this project

**cd** - change directory (navigate to the place where you want to create the file)

**mnt** - mount (we work in Windows but use Linux; this helps navigate to the place where we want to go. It makes the location available within the Linux file system.)
```bash
cd /mnt/c/Users/Kejty
```
* note: if in your pathway to the desired location are spaces in names, e.g. /mnt/c/Users/Kejty laptop! Backslash (\) needs to be inserted instead of the space, with one space after the slash
and the location needs to be wrapped in double quote = "/mnt/c/Users/Kejty\ laptop"        -to avoid it alsways use underscore _

**mkdir** - make directory (the file we want to create)
```bash
mkdir -p Tutorials        #- p stands for parent; we can create multiple folders at once. Kejty is the parent folder. If it was missing, we could create mkdir -p Kejty/Tutorials - it would create the parent folder (Kejty), and inside the folder, the folder Tutorials
```
We have our metadata -> if not, please 

Now, creating a manifest.

Manifest tells Qiime2 which FASTQ files belong to which sample.

🗺️ Manifest is like a map for Qiime 2 between the sample names and two FASTQ files (forward and reverse)

So, for example:

```text
sample-id    forward-file              reverse-file
sample1      sample1_R1.fastq.gz       sample1_R2.fastq.gz
sample2      sample2_R1.fastq.gz       sample2_R2.fastq.gz
```
 **fastq_dir=** will create a path to the directory/folder where the FASTQ files are saved.
```bash
fastq_dir=/mnt/c/Users/Kejty
```
 
