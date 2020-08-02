set editconf "gmx editconf" 
set grompp "gmx grompp" 
set make_ndx "gmx make_ndx" 
# FILES INPUTS 
set filename ionized 
# Write Gromacs Compatible Top and tpr 
# WARNING JUST FOR ANALISIS !!!!!!! 
if { [file exists ${filename}.top] == 0} { 
package require topotools 
topo writegmxtop ${filename}.top 
} 
## Gromacs Files 
# gro & tpr 
if { [file exists ${filename}.gro] == 0 || [file exists ${filename}.tpr] == 0} { 
            catch {exec $editconf -f ${filename}.pdb -o ${filename}.gro -nopbc} 
            catch {exec $grompp -f ions.mdp -c ${filename}.gro -maxwarn 10 -p ${filename}.top -o ${filename}.tpr } 
    }
