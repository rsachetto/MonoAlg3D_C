# MonoAlg3D ![build](https://github.com/rsachetto/MonoAlg3D_C/actions/workflows/build.yml/badge.svg)

The MonoAlg3D is a program for solving the 3D monodomain equation by applying the Finite Volume Method.

# Pre-Requisites

  - Linux-based operating system
  - Nvidia Driver 
  - CUDA library

### Setting the enviroment

Ubuntu: Refer to [the ubuntu guide](guide-monoalg3d-ubuntu.md)

Fedora: Refer to [the fedora guide](guide-monoalg3d-fedora.md)

Windows: Refer to [the windows guide](guide-monoalg3d-windows.md)

### Compile
```sh
$ ./build.sh
```
The binary files will be saved in the ```bin``` folder.

# Running examples
```sh
$ bin/MonoAlg3D -c example_configs/cuboid_ohara.ini 
```

The output will be saved in the VTK format. In order to see the results you can use Paraview (https://www.paraview.org/) or the compiled visualization tool in ```bin/MonoAlg3D_visualizer```. You can also set the output to plain text, by changing the section ```save_result``` in example_configs/cuboid_ohara.ini to:

```ìni
[save_result]
print_rate=250
output_dir=./outputs/tmp_cube
main_function=save_as_text_or_binary
binary=false
```

In the plain text format we have:

- Each line represents a Volume
- Each volume is represented by its center point (X, Y, and Z), the value of half of its side length on x, y and z and the calculated V

Example file:

```
850,850,950,50,50,50, -85
850,950,950,50,50,50, -85
850,950,850,50,50,50, -85
```

This file represents 3 volumes with 100 micrometer of side. The first volume is centered at  at 850,850,950 and the calculated V is -85 mV.

# Contributors:

@rsachetto Rafael Sachetto Oliveira

@bergolho Lucas Arantes Berg

@Rodrigo-Weber-dos-Santos Rodrigo Weber dos Santos

Among others.

# How to cite:

Oliveira RS, Rocha BM, Burgarelli D, Meira Jr W, Constantinides C, dos Santos RW. Performance evaluation of GPU parallelization, space‐time adaptive algorithms, and their combination for simulating cardiac electrophysiology. Int J Numer Meth Biomed Engng. 2018;34:e2913. https://doi.org/10.1002/cnm.2913

# Credits
[Heart icons created by phatplus - Flaticon](https://www.flaticon.com/free-icons/heart)

Of course. Here is the appendix formatted as a Markdown file with the title, authors, and a link for each publication.

# Research utilizing MonoAlg3D (not all, not sorted)

### In Silico TRials guide optimal stratification of ATrIal FIbrillation patients to Catheter Ablation and pharmacological medicaTION: the i-STRATIFICATION study

  * **Authors:** Albert Dasí, Claudia Nagel, Michael T B Pope, Rohan S Wijesurendra, Timothy R Betts, Rafael Sachetto, Axel Loewe, Alfonso Bueno-Orovio, Blanca Rodriguez
  * **Link:** [https://doi.org/10.1093/europace/euae150](https://doi.org/10.1093/europace/euae150)

### Mechanisms of ischaemia-induced arrhythmias in hypertrophic cardiomyopathy: a large-scale computational study

  * **Authors:** James A Coleman, Ruben Doste, Zakariye Ashkir, Raffaele Coppini, Rafael Sachetto, Hugh Watkins, Betty Raman, Alfonso Bueno-Orovio
  * **Link:** [https://doi.org/10.1093/cvr/cvae086](https://doi.org/10.1093/cvr/cvae086)

### Uncovering Arrhythmic Substrate in a Heart Failure Patient Using Digital Twins

  * **Authors:** Thaís de J. Soares, Guilherme M. Couto, Yan B. Werneck, Daniel K. Almeida, Daniel M. P. Leme, João P. B. Pereira, Filipe L. Namorato, Matheus C. Faesy, Tiago D. Franco, Fabricio J. M. Santos, Raul P. Barra, Marcelle C. S. B. Vasconcelos, Thaiz R. Schmal, Thiago G. S. Souza, Rafael S. Oliveira, Joventino O. Campos, and Rodrigo W. dos Santos
  * **Link:** [https://cinc.org/2025/Program/accepted/301\_Preprint.pdf](https://cinc.org/2025/Program/accepted/301_Preprint.pdf)

### Variability in electrophysiological properties and conducting obstacles controls re-entry risk in heterogeneous ischaemic tissue

  * **Authors:** Brodie A. J. Lawson, Rafael S. Oliveira, Lucas A. Berg, Pedro A. A. Silva, Kevin Burrage, Rodrigo Weber dos Santos
  * **Link:** [https://doi.org/10.1098/rsta.2019.0341](https://www.google.com/search?q=https://doi.org/10.1098/rsta.2019.0341)

### Ectopic beats arise from micro-reentries near infarct regions in simulations of a patient-specific heart model

  * **Authors:** Rafael Sachetto Oliveira, Sergio Alonso, Fernando Otaviano Campos, Bernardo Martins Rocha, João Filipe Fernandes, Titus Kuehne, Rodrigo Weber Dos Santos
  * **Link:** [https://doi.org/10.1038/s41598-018-34304-y](https://doi.org/10.1038/s41598-018-34304-y)

### Killing many birds with two stones: hypoxia and fibrosis can generate ectopic beats in a human ventricular model

  * **Authors:** Rafael Sachetto, Sergio Alonso, Rodrigo Weber dos Santos
  * **Link:** [https://doi.org/10.3389/fphys.2018.00764](https://doi.org/10.3389/fphys.2018.00764)

### Digital twinning of the human ventricular activation sequence to Clinical 12-lead ECGs and magnetic resonance imaging using realistic Purkinje networks for in silico clinical trials

  * **Authors:** J. Camps, L.A. Berg, Z.J. Wang, R. Sebastian, L.L. Riebel, R. Doste, X. Zhou, et al.
  * **Link:** [https://doi.org/10.1016/j.media.2024.103108](https://doi.org/10.1016/j.media.2024.103108)

### Development of a Three-Dimensional Computational Pipeline in Python for Personalized Heart Modeling

  * **Authors:** Filipe L Namorato, Daniel M P Leme, Thaıs J Soares, Rafael S Oliveira, Thaiz R Schmal, Rodrigo W dos Santos, Joventino O Campos
  * **Link:** [https://cinc.org/2025/Program/accepted/329\_Preprint.pdf](https://cinc.org/2025/Program/accepted/329_Preprint.pdf)

### In-silico drug trials for precision medicine in atrial fibrillation: From ionic mechanisms to electrocardiogram-based predictions in structurally-healthy human atria

  * **Authors:** A. Dasí, A. Roy, R. Sachetto, J. Camps, A. Bueno-Orovio, B. Rodriguez
  * **Link:** [https://doi.org/10.3389/fphys.2022.966046](https://doi.org/10.3389/fphys.2022.966046)
  
### Benchmarking Open Cardiac Electrophysiology Simulators: MonoAlg3D and OpenCARP

  * **Authors:** de Lima, Lucas MR and Oliveira, Rafael S and Campos, Fernando O and Berg, Lucas A and Campos, Joventino O and dos Santos, Rodrigo W
  * **Link:** [https://cinc.org/2025/Program/accepted/376.html](https://cinc.org/2025/Program/accepted/376.html) [17]
 
### ECG analysis of ventricular fibrillation dynamics reflects ischaemic progression subject to variability in patient anatomy and electrode location

  * **Authors:** Hector Martinez-Navarro, Ambre Bertrand, Ruben Doste, Hannah Smith, Jakub Tomek, Giuseppe Ristagno, Rafael S. Oliveira, Rodrigo Weber dos Santos, and Blanca Rodriguez
  * **Link:** [https://www.frontiersin.org/journals/cardiovascular-medicine/articles/10.3389/fcvm.2024.1408822/full](https://www.frontiersin.org/journals/cardiovascular-medicine/articles/10.3389/fcvm.2024.1408822/full)

### Creating Digital Twins for Cardiac Electrophysiology in Patients with Non-Ischemic Dilated Cardiomyopathy

  * **Authors:** Thais de Jesus Soares, Tiago Dutra Franco, Fabricio Junio Mendes Santos, Thaiz Ruberti Schmal, Thiago Goncalves Schroder e Souza, Rafael Sachetto Oliveira, Joventino de Oliveira Campos, Rodrigo Weber dos Santos
  * **Link:** [https://doi.org/10.1007/978-3-031-94921-0\_64](https://www.google.com/search?q=https://doi.org/10.1007/978-3-031-94921-0_64)
  

### Electrophysiological mechanisms underlying T wave pseudonormalisation on stress ECGs in hypertrophic cardiomyopathy

  * **Authors:** J. Coleman, R. Doste, M. Beltrami, I. Olivotto, B. Raman, A. Bueno Orovio
  * **Link:** [https://doi.org/10.1016/j.compbiomed.2023.107829](https://www.google.com/search?q=https://doi.org/10.1016/j.compbiomed.2023.107829)
  

### A Multi-Scale Computational Framework for Human-Based Modelling and Simulation of Adverse Cardiac Remodelling in Type 2 Diabetes

  * **Authors:** Ambre Bertrand, Lucas Arantes-Berg, Ruben Doste, Jakub Tomek, Albert Dasi, Abhirup Banerjee, Julia Camps, Vicente Grau, Blanca Rodriguez
  * **Link:** [https://cinc.org/2025/Program/accepted/418.pdf](https://cinc.org/2025/Program/accepted/418.pdf)

