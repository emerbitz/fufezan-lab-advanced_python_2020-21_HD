import pandas as pd


def read_msp_file(file_path='../data/cptac2_mouse_hcd_selected.msp',
                  mz_min=300,
                  mz_max=900,
                  collision_energy_min=25,
                  collision_energy_max=35,
                  sequence_length_max=30
                  ):
    """Extracts the spectra and sequences out of a msp-file.

    Args:
        file_path: Relative path to msp-file as string
        mz_min: Minimal value for m/z.
        mz_max: Maximal value for m/z.
        collision_energy_min: Minimal value for the collision energy.
        collision_energy_max: Maximal value for the collision energy.
        sequence_length_max: Maximal length of the amino acide sequence.

    Returns: Spectra as dataframe, sequence as dataframe

    """

    spectrum_df = pd.DataFrame()
    sequence_lst = []
    with open(file_path) as msp_file:
        for box_num, box in enumerate(msp_file.read().split('\n\n')[:-1]):
            print(box_num)
            head, tail = box.split('\n', maxsplit=1)
            seq, col_energy = head[6:-2].split('/')
            col_energy = float(col_energy.split('_')[-1])
            #         print(head)
            #         print(seq)
            #         print(col_energy)
            #         print(len(seq))
            if ((collision_energy_min <= col_energy <= collision_energy_max) & (
                    len(seq) <= sequence_length_max)) is False:
                continue
            #         print(box_num)
            sequence_lst.append(seq)
            spectrum_lst = []
            for pos, line in enumerate(tail.split('\n')[3:]):
                #             print(pos)
                #             print(line)
                mz, intensity = line.split('\t')[:2]
                mz = round(float(mz))
                if (mz_min <= mz <= mz_max) is False:
                    continue
                intensity = float(intensity)
                spectrum_lst.append((mz, intensity))
            spectrum_dic = dict(sorted(spectrum_lst))
            spectrum_df = spectrum_df.append(spectrum_dic, ignore_index=True)
            # if box_num == 5:
            #     break
    sequence_df = pd.DataFrame(sequence_lst, columns=['sequence'])
    mask = ~(spectrum_df == 0).all(axis=1)
    spectrum_df = spectrum_df[mask]
    return spectrum_df, sequence_df


if __name__ == '__main__':
    spectrum, sequence = read_msp_file('../data/cptac2_mouse_hcd_selected.msp')
    print(sequence)
    print(spectrum)
