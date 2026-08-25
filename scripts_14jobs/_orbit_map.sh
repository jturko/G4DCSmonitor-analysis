# _orbit_map.sh -- sourced by generate_run_macros.sh and submit_all.sh.
#
# Forward orbit map for the CASTOR 440/84 basket under cask rotation about z.
# Convention: /dcs-monitor/det/setRotation 0 0 <deg> applies an active CCW
# rotation about +z, so e.g. base fuel 0 + 60 deg -> global fuel 28.
#
# Provides two functions:
#   lookup_global_fuel <base_fuel> <rotation_deg>   -> echoes global fuel
#   lookup_base_rot    <global_fuel>                -> echoes "BASE ROT"
#
# Each integer 0..83 appears exactly once across (base, rot).

lookup_global_fuel() {
    local base="$1" rot="$2"
    case "${base}:${rot}" in
        0:0)    echo  0 ;;   0:60)   echo 28 ;;   0:120)  echo 73 ;;
        0:180)  echo 83 ;;   0:240)  echo 55 ;;   0:300)  echo 10 ;;

        1:0)    echo  1 ;;   1:60)   echo 19 ;;   1:120)  echo 65 ;;
        1:180)  echo 82 ;;   1:240)  echo 64 ;;   1:300)  echo 18 ;;

        2:0)    echo  2 ;;   2:60)   echo 11 ;;   2:120)  echo 56 ;;
        2:180)  echo 81 ;;   2:240)  echo 72 ;;   2:300)  echo 27 ;;

        3:0)    echo  3 ;;   3:60)   echo  4 ;;   3:120)  echo 46 ;;
        3:180)  echo 80 ;;   3:240)  echo 79 ;;   3:300)  echo 37 ;;

        5:0)    echo  5 ;;   5:60)   echo 38 ;;   5:120)  echo 74 ;;
        5:180)  echo 78 ;;   5:240)  echo 45 ;;   5:300)  echo  9 ;;

        6:0)    echo  6 ;;   6:60)   echo 29 ;;   6:120)  echo 66 ;;
        6:180)  echo 77 ;;   6:240)  echo 54 ;;   6:300)  echo 17 ;;

        7:0)    echo  7 ;;   7:60)   echo 20 ;;   7:120)  echo 57 ;;
        7:180)  echo 76 ;;   7:240)  echo 63 ;;   7:300)  echo 26 ;;

        8:0)    echo  8 ;;   8:60)   echo 12 ;;   8:120)  echo 47 ;;
        8:180)  echo 75 ;;   8:240)  echo 71 ;;   8:300)  echo 36 ;;

        13:0)   echo 13 ;;   13:60)  echo 39 ;;   13:120) echo 67 ;;
        13:180) echo 70 ;;   13:240) echo 44 ;;   13:300) echo 16 ;;

        14:0)   echo 14 ;;   14:60)  echo 30 ;;   14:120) echo 58 ;;
        14:180) echo 69 ;;   14:240) echo 53 ;;   14:300) echo 25 ;;

        15:0)   echo 15 ;;   15:60)  echo 21 ;;   15:120) echo 48 ;;
        15:180) echo 68 ;;   15:240) echo 62 ;;   15:300) echo 35 ;;

        22:0)   echo 22 ;;   22:60)  echo 40 ;;   22:120) echo 59 ;;
        22:180) echo 61 ;;   22:240) echo 43 ;;   22:300) echo 24 ;;

        23:0)   echo 23 ;;   23:60)  echo 31 ;;   23:120) echo 49 ;;
        23:180) echo 60 ;;   23:240) echo 52 ;;   23:300) echo 34 ;;

        32:0)   echo 32 ;;   32:60)  echo 41 ;;   32:120) echo 50 ;;
        32:180) echo 51 ;;   32:240) echo 42 ;;   32:300) echo 33 ;;

        *)
            echo "[orbit_map] no entry for (base=$base, rot=$rot)" >&2
            return 2 ;;
    esac
}

# Inverse: global fuel index in 0..83 -> "BASE ROT"
lookup_base_rot() {
    local g="$1"
    case "$g" in
         0) echo  "0   0" ;;    1) echo  "1   0" ;;    2) echo  "2   0" ;;
         3) echo  "3   0" ;;    4) echo  "3  60" ;;    5) echo  "5   0" ;;
         6) echo  "6   0" ;;    7) echo  "7   0" ;;    8) echo  "8   0" ;;
         9) echo  "5 300" ;;   10) echo  "0 300" ;;   11) echo  "2  60" ;;
        12) echo  "8  60" ;;   13) echo "13   0" ;;   14) echo "14   0" ;;
        15) echo "15   0" ;;   16) echo "13 300" ;;   17) echo  "6 300" ;;
        18) echo  "1 300" ;;   19) echo  "1  60" ;;   20) echo  "7  60" ;;
        21) echo "15  60" ;;   22) echo "22   0" ;;   23) echo "23   0" ;;
        24) echo "22 300" ;;   25) echo "14 300" ;;   26) echo  "7 300" ;;
        27) echo  "2 300" ;;   28) echo  "0  60" ;;   29) echo  "6  60" ;;
        30) echo "14  60" ;;   31) echo "23  60" ;;   32) echo "32   0" ;;
        33) echo "32 300" ;;   34) echo "23 300" ;;   35) echo "15 300" ;;
        36) echo  "8 300" ;;   37) echo  "3 300" ;;   38) echo  "5  60" ;;
        39) echo "13  60" ;;   40) echo "22  60" ;;   41) echo "32  60" ;;
        42) echo "32 240" ;;   43) echo "22 240" ;;   44) echo "13 240" ;;
        45) echo  "5 240" ;;   46) echo  "3 120" ;;   47) echo  "8 120" ;;
        48) echo "15 120" ;;   49) echo "23 120" ;;   50) echo "32 120" ;;
        51) echo "32 180" ;;   52) echo "23 240" ;;   53) echo "14 240" ;;
        54) echo  "6 240" ;;   55) echo  "0 240" ;;   56) echo  "2 120" ;;
        57) echo  "7 120" ;;   58) echo "14 120" ;;   59) echo "22 120" ;;
        60) echo "23 180" ;;   61) echo "22 180" ;;   62) echo "15 240" ;;
        63) echo  "7 240" ;;   64) echo  "1 240" ;;   65) echo  "1 120" ;;
        66) echo  "6 120" ;;   67) echo "13 120" ;;   68) echo "15 180" ;;
        69) echo "14 180" ;;   70) echo "13 180" ;;   71) echo  "8 240" ;;
        72) echo  "2 240" ;;   73) echo  "0 120" ;;   74) echo  "5 120" ;;
        75) echo  "8 180" ;;   76) echo  "7 180" ;;   77) echo  "6 180" ;;
        78) echo  "5 180" ;;   79) echo  "3 240" ;;   80) echo  "3 180" ;;
        81) echo  "2 180" ;;   82) echo  "1 180" ;;   83) echo  "0 180" ;;
        *) echo "[orbit_map] global fuel out of range: $g" >&2; return 2 ;;
    esac
}

