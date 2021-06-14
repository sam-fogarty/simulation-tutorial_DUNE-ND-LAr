import argparse
import time
import json

import larpix

_default_vdda = 1800
_default_vref_dac = 217
_default_vcm_dac = 71

def key2unique(key, channel):
    io_group,io_channel,chip_id = str(key).split('-')
    return ((int(io_group)*256 + int(io_channel))*256 + int(chip_id))*64 + int(channel)

def dac2mv(dac, max, bits=8):
    return max * dac/(2**bits)

def main(controller_config, vdda=_default_vdda, vref_dac=_default_vref_dac,
    vcm_dac=_default_vcm_dac, **kwargs):
    c = larpix.Controller()
    c.load(controller_config)

    config_dict = dict()
    for chip_key in c.chips:
        for channel in range(64):
            unique = key2unique(chip_key,channel)
            vref_mv = dac2mv(vref_dac,vdda)
            vcm_mv = dac2mv(vcm_dac,vdda)
            config_dict[unique] = dict(
                vref_mv=vref_mv,
                vcm_mv=vcm_mv,
                )
    dt = time.strftime('%y-%m-%d_%H-%M-%S')
    with open('evd_config_{}.json'.format(dt),'w') as fo:
        json.dump(config_dict, fo, sort_keys=True, indent=4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--controller_config','-c',required=True,type=str)
    parser.add_argument('--vdda',default=_default_vdda,type=float,help='''default=%(default)s mV''')
    parser.add_argument('--vref_dac',default=_default_vref_dac,type=int,help='''default=%(default)s''')
    parser.add_argument('--vcm_dac',default=_default_vcm_dac,type=int,help='''default=%(default)s''')
    args = parser.parse_args()
    main(**vars(args))
