summercamp2014
==============

The ntuples are on the LEPP `lnx` cluster:

```
[jmt46@lnx201 ~] $ ls -lh /cda/cdat/tem/jmt46/flattrees/
total 25G
drwxr-xr-x 2 jmt46 cms 4.0K Jun 17 10:11 logs
-rw-r--r-- 1 jmt46 cms 196M Jun 17 01:03 mfv_neutralino_tau0100um_M0400.root
-rw-r--r-- 1 jmt46 cms 213M Jun 17 01:04 mfv_neutralino_tau0300um_M0400.root
-rw-r--r-- 1 jmt46 cms 231M Jun 17 01:04 mfv_neutralino_tau1000um_M0400.root
-rw-r--r-- 1 jmt46 cms 244M Jun 17 01:05 mfv_neutralino_tau9900um_M0400.root
-rw-r--r-- 1 jmt46 cms  60M Jun 17 01:05 qcdht0100.root
-rw-r--r-- 1 jmt46 cms 1.2G Jun 17 01:06 qcdht0250.root
-rw-r--r-- 1 jmt46 cms 6.9G Jun 17 01:08 qcdht0500.root
-rw-r--r-- 1 jmt46 cms 4.6G Jun 17 01:11 qcdht1000.root
-rw-r--r-- 1 jmt46 cms 1.3G Jun 17 01:12 ttbardilep.root
-rw-r--r-- 1 jmt46 cms 4.0G Jun 17 01:16 ttbarhadronic.root
-rw-r--r-- 1 jmt46 cms 5.5G Jun 17 01:19 ttbarsemilep.root
```

You can copy them to your laptop, but be prepared for it to take a
while. The syntax to copy is

```
  scp username@lnx201.lns.cornell.edu:/cdat/tem/jmt46/ttbarhadronic.root .
```

Some skeleton code to read them and output the number of events
consists of the `.C` and `.h` file in this directory.

Put both the `.C` and `.h` in the same directory as the `.root` files,
and then you can do

```
[user@localhost $]  root -l MFVFlatNtupleReader.C+
root [0] MFVFlatNtupleReader r
root [1] r.Loop()
```

Be sure not to forget the `+` so that root compiles the macro, or it
will be slower than it already is. You should eventually see output
like

```
done!
sample mfv_neutralino_tau0100um_M0400 events read:      96043
sample mfv_neutralino_tau0300um_M0400 events read:      96653
sample mfv_neutralino_tau1000um_M0400 events read:      96975
sample mfv_neutralino_tau9900um_M0400 events read:      96734
sample                  ttbarhadronic events read:    2503538
sample                   ttbarsemilep events read:    3489937
sample                     ttbardilep events read:     796885
sample                      qcdht0100 events read:      37489
sample                      qcdht0250 events read:     822457
sample                      qcdht0500 events read:    4759413
sample                      qcdht1000 events read:    3088150
```
