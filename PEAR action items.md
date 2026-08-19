# GitHub repository

| To do |  |  |  |
| :---- | :---- | :---- | :---- |
| ![No type][image1] Task | ![Drop-downs][image2] Status | ![Drop-downs][image2]Priority | ![No type][image1] Notes |
| Ensure correct package versions installed | Launched | High | NumPy \< 2.0 required |
| Document versions used | Launched | Medium | Run 'pip freeze \> versions.txt' or similar |
| Improve error messages shown when hashrf fails | Launched | High | Print the errors thrown by the executable, which would be helpful for debugging in the following cases: hashrf failing due to lack of permissions (would print “Permission denied”). Incorrect number of trees being passed to hashrf executable due to missing newline characters (would print “\*\*\* Number of trees in the input file: 1 Fatal error: at least two trees expected.”) |
| Fix issue with concatenated Newick files | Launched | High | Use a more robust method of counting the number of trees in a file, or force a newline to be inserted following each Newick string. When trees with no linefeed at the end of the string are passed to PEAR, all trees in the Set\_collection end up in a single line, which leads to the wrong number of trees being passed to the hashrf executable. Opening the file in an editor that appends an \\x0a to end upon saving it appears to fix the problem. Ensuring that the code that processes the Newick files enforces this ought to fix the issue. |
| Ensure hashrf has correct permissions and is executable | Launched | Lower | Notebook examples did not run, until the permissions were adjusted (possibly because running from within examples folder leads to local files being imported) |
| Ensure no anaconda channels used | Launched | Medium | EBI no longer allows the use of Anaconda’s default channels (this is because research institutes are treated differently than universities by Anaconda Inc.). Make sure that everything runs using packages from conda-forge. Miniforge should do this by default. |

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAQAQMAAAAs1s1YAAAABlBMVEUAAABER0byc6G0AAAAAXRSTlMAQObYZgAAAB9JREFUeF5jYEAD9h8YmEA0MwOYZmSWWQjhs4H56BgAT4ECDeGaeV4AAAAASUVORK5CYII=>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAQCAYAAAAWGF8bAAAAx0lEQVR4Xu2TYRHCMAyFKwEJSEBCjyVpXIAEHIATJCBhEpCAhEkA0tEtTVcod/zku8ufvDR7fduc+/NTmHkNge6t1QU82x0TQHRMg8i4t7rGI266QJc0b3Xt7Gq1d3jvV69zQyZUn9QAMHg5K8vnpuTBtFNzXzEawj5rKL0AmE7pFku3KXrR8jPHeaRELy00m2NhuYIstT1hPB8OqoF9dKmDbQQD3ZZcTznIJ2S1GmlZ5k4DhENa3FrbDz+Bk5eTIqhVdFbJ8wG0lJX5M/zhmwAAAABJRU5ErkJggg==>