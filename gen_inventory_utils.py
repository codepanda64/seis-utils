from obspy import read_inventory, Inventory, UTCDateTime
from obspy.core.inventory import Network, Station, Channel, Site, Response
from obspy import read_inventory
from obspy.clients.nrl import NRL

from datetime import datetime
from enum import Enum


class InstrumentType(Enum):
    """目前支持的仪器类型"""
    
    REFTEK130B_40T2S = 1, 'REFTEK130B+CMG-40T2S_100'
    REFTEK130B_40T30S = 2, 'REFTEK130B+CMG-40T30S_100'
    REFTEK130S_40T2S = 3, 'REFTEK130S+CMG-40T2S_100'
    REFTEK130S_40T30S = 4, 'REFTEK130S+CMG-40T30S_100'
    REFTEK130B_40T2S_200 = 5, 'REFTEK130B+CMG-40T2S_200'

class GenerateInventoryUtil():
    
    def __init__(self, stationxml: str=None, default_end_date:str='2599-01-01T00:00:00.000000Z') -> None:
        if stationxml:
            self.inv = read_inventory(stationxml)
        else:
            self.inv = Inventory()
        self.default_end_date = UTCDateTime(default_end_date)
        self.nrl = NRL()
    
    
    def _response_generator(self, instrument_keys: InstrumentType) -> Response:
        """ 根据关键字获取仪器相应

        Args:
            instrument_keys (InstrumentType): 仪器类型关键字

        Returns:
            Response: 仪器相应
        """
        if instrument_keys == InstrumentType.REFTEK130B_40T2S:
            response = self.nrl.get_response(
                sensor_keys=['Guralp', 'CMG-40T', '2s - 100Hz', '2000'],
                datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '100']
            )
        elif instrument_keys == InstrumentType.REFTEK130B_40T30S:
            response = self.nrl.get_response(
                sensor_keys=['Guralp', 'CMG-40T', '30s - 100Hz', '2000'],
                datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '100']
            )
        elif instrument_keys == InstrumentType.REFTEK130S_40T2S:
            response = self.nrl.get_response(
                sensor_keys=['Guralp', 'CMG-40T', '2s - 100Hz', '2000'],
                datalogger_keys=['REF TEK', 'RT 130S & 130-SMHR', '1', '100']
            )
        elif instrument_keys == InstrumentType.REFTEK130S_40T30S:
            response = self.nrl.get_response(
                sensor_keys=['Guralp', 'CMG-40T', '30s - 100Hz', '2000'],
                datalogger_keys=['REF TEK', 'RT 130S & 130-SMHR', '1', '100']
            )
        elif instrument_keys == InstrumentType.REFTEK130B_40T2S_200:
            response = self.nrl.get_response(
                sensor_keys=['Guralp', 'CMG-40T', '2s - 100Hz', '2000'],
                datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200']
            )
            
        return response
    
    
    def change_latest_location(self, network_code: str, station_code: str, latitude: float, longitude: float, elevation: float, changed: UTCDateTime, depth: float=0.0) -> None:
        """台站位置信息发生更改时，增加新的位置到inventory文件中

        Args:
            network_code (str): 台网代码
            station_code (str): 台站代码
            latitude (float): 台站新的纬度
            longitude (float): 台站新的经度
            elevation (float): 台站新的高程
            changed (UTCDateTime): 位置变更时间，北京时间， 例如：2023-03-31T17:41:00+08'
            depth (float): 台站新的深度
            

        Raises:
            Exception: 台站在inventory中找不到时抛出异常
        """
        latest_start_date = changed
        previous_end_date = latest_start_date - 2 * 3600
        
        for i, network in enumerate(self.inv.networks):
            if network.code == network_code:
                for j, station in enumerate(self.inv.networks[i].stations):
                    if station.code == station_code and (station.end_date == self.default_end_date or station.end_date is None):
                        new_station = station.copy()
                        self.inv.networks[i].stations[j].end_date = previous_end_date
                        
                        new_station.start_date = latest_start_date
                        new_station.latitude = latitude
                        new_station.longitude = longitude
                        new_station.elevation = elevation
                        for ch in new_station.channels:
                            ch.latitude = latitude
                            ch.longitude = longitude
                            ch.elevation = elevation
                        self.inv.networks[i].stations.append(new_station)
                        break
        else:
            raise Exception("The station was not found in this inventory!!!")
        
    
    def change_latest_instrument(self, network_code: str, station_code: str, instrument_keys: InstrumentType, changed: UTCDateTime) -> None:
        """当台站仪器发生变更时更新inventory文件

        Args:
            network_code (str): 台网代码
            station_code (str): 台站代码
            instrument_keys (InstrumentType): 仪器关键字符串（用来获取仪器相应）
            changed (UTCDateTime): 变更时间， 北京时间， 例如：2023-03-31T17:41:00+08'
            
        Raises:
            Exception: 台站在inventory中找不到时抛出异常
        """
        latest_start_date = changed
        previous_end_date = latest_start_date - 2*3600
        
        for i, network in enumerate(self.inv.networks):
            if network.code == network_code:
                for j, station in enumerate(self.inv.networks[i].stations):
                    if station.code == station_code and (station.end_date == self.default_end_date or station.end_date is None):
                        new_station = station.copy()
                        self.inv.networks[i].stations[j].end_date = previous_end_date
                        
                        new_station.start_date = latest_start_date
                        for ch in new_station.channels:
                            ch.response = self._response_generator(instrument_keys=instrument_keys)
                            ch.response = self._response_generator(instrument_keys=instrument_keys)
                            ch.response = self._response_generator(instrument_keys=instrument_keys)
                            
                        self.inv.networks[i].stations.append(new_station)
                        break
                    
        else:
            raise Exception("The station was not found in this inventory!!!")        
        
    
    def create_station_or_change_latest_station_all(self, network_code: str, station_code: str, latitude: float, longitude: float, elevation: float, instrument_keys: InstrumentType, created: UTCDateTime, depth: float=0.0, sample_rate: int=100, name: str=None, is_create: bool=True) -> None:
        """新添加台站或者台站搬迁并且更换仪器

        Args:
            network_code (str): 台网代码
            station_code (str): 台站代码
            latitude (float): 台站纬度
            longitude (float): 台站经度
            elevation (float): 台站高程
            instrument_keys (InstrumentType): 台站仪器关键字
            created (UTCDateTime): 创建日期，北京时间， 例如：2023-03-31T17:41:00+08'
            depth (float, optional): 台站埋深. Defaults to 0.0.
            sample_rate (int, optional): 采样率. Defaults to 100.
            name (str, optional): 台站名. Defaults to None.
            is_create (bool, optional): 是否是创建新台站，如果为否则视为台站搬迁并更换仪器. Defaults to True.

        Raises:
            Exception: 如果台站已存在抛出异常
        """
        
        start_date = created
        new_station = Station(code=station_code, latitude=latitude, longitude=longitude, elevation=elevation, creation_date=start_date, start_date=start_date, end_date=self.default_end_date, site=Site(name=name))
        response = self._response_generator(instrument_keys=instrument_keys)
        for ch in ['SHZ', 'SHN', 'SHE']:
            if ch == 'SHZ':
                az = 0.0
                dip = -90.0
            elif ch == 'SHN':
                az = 0.0
                dip = 0.0
            elif ch == 'SHE':
                az = 90.0
                dip = 0.0
            
            new_channel = Channel(code=ch, location_code='00', latitude=latitude, longitude=longitude, elevation=elevation, depth=depth, azimuth=az, dip=dip, sample_rate=sample_rate)
            new_channel.response = response
            new_station.channels.append(new_channel)

        for i, network in enumerate(self.inv.networks):
            if network.code == network_code:
                for j, station in enumerate(self.inv.networks[i].stations):
                    if station.code == station_code:
                        if is_create:
                            raise Exception("The station already exists!!!")
                        else:
                            previous_end_date = start_date - 2*3600
                            self.inv.networks[i].stations[j].end_date = previous_end_date
                            self.inv.networks[i].stations.append(new_station)
                            break
                else:
                    self.inv.networks[i].stations.append(new_station)
                    break
        else:
            new_network = Network(code=network_code)
            new_network.stations.append(new_station)
            self.inv.networks.append(new_network)
        
            
    def write(self, filepath=None):
        if not filepath:
            filepath = f'station.{datetime.now().strftime("%Y%m%d%H%M%S")}.xml'
        self.inv.write(filepath, format="stationxml", validate=True)
        

if __name__ == '__main__':
    import sys
    import argparse
    from argparse import RawTextHelpFormatter
    
    instruments_info_str = ""
    for instr in InstrumentType:
        instruments_info_str += f'\n {instr.value[0]}: {instr.value[1]}'
    
    parser = argparse.ArgumentParser(description='Generate Inventory Util', formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('--station-code', type=str,  help='station code. format: NET.STA, e.g:G1.53251')
    parser.add_argument('--station-name', type=str,  help='station name')
    
    parser.add_argument('--is-create', action='store_true', default=False, help='create a new station into the inventory')
    
    parser.add_argument('--is-change-loc', action='store_true', default=False, help='modifiy a station location')
    parser.add_argument('--lat', type=float, help='The station new latitude')
    parser.add_argument('--lon', type=float, help='The station new longitude')
    parser.add_argument('--elev', type=float, help='The station new elevation')
    parser.add_argument('--dep', type=float, default=0.0, help='The station new depth')
    
    parser.add_argument('--is-change-ins', action='store_true', default=False, help='Changing the station\'s instrument')
    parser.add_argument('--ins-key', type=int, help=f'Instrument keywords for station: {instruments_info_str}')
    
    parser.add_argument('--change-date', type=str, help='Change date. format: YYYY-mm-ddTHH:MMZ, e.g: 2023-05-10T10:20+08')
    parser.add_argument('--station-xml', type=str, help='Station XML file')
    
    args = parser.parse_args()
        
    # station_code = args.station_code
    if args.station_code is not None:
        try:
            station_code = args.station_code.split('.')
            net_code = station_code[0].upper()
            sta_code = station_code[1].upper()
        except Exception as e:
            print(f'--station-code "{args.station_code}" is Wrong format. format: NET.STA, e.g: G1.53251. {e}')
    station_name = args.station_name

    is_create = args.is_create
    
    is_change_loc = args.is_change_loc
    is_change_ins = args.is_change_ins
    lat = args.lat
    lon = args.lon
    elev = args.elev
    dep = args.dep
    ins_key = args.ins_key
    stationxml = args.station_xml
    change_date = None
    if args.change_date is not None:
        try:
            change_date = UTCDateTime(args.change_date)
        except Exception as e:
            print(f"--change-date: date format is incorrect. format: YYYY-mm-ddTHH:MMZ, e.g. 2023-05-10T08:40+08 \n{e}")
    
    if not (is_create or is_change_loc or is_change_ins):
        parser.print_help()
        sys.exit()
    
    
    if is_create and (is_change_loc or is_change_ins):
        print("################################################################################################################################\n" +
              "####  !!!Creating a new station cannot be True together with changing the location and instrumentation of an old station!!! ####\n" +
              "################################################################################################################################\n")
        parser.print_help()
        sys.exit()
    
    inventory_gen = GenerateInventoryUtil(stationxml=stationxml)
    
    if is_create or (is_change_loc and is_change_ins):
        if sta_code is None or lat is None or lon is None or elev is None or ins_key is None or change_date is None:
            print("#######################################################################\n" +
                  "#### !!!Creating a new station must set lat, lon, elev, ins_key!!! ####\n" +
                  "#######################################################################\n")
            parser.print_help()
            sys.exit()
        else:
            for instrument in InstrumentType:
                if ins_key == instrument.value[0]:
                    print(instrument)
                    instrument_str = instrument.value[1]
                    inventory_gen.create_station_or_change_latest_station_all(network_code=net_code, station_code=sta_code, latitude=lat, longitude=lon, elevation=elev, instrument_keys=instrument, created=change_date, name=station_name, is_create=is_create, depth=dep)
                    break
            else:
                print("##########################################################\n" +
                      f"#### !!!ins-key must set {instruments_info_str} !!! #####\n" +
                      "##########################################################\n")
                parser.print_help()
                sys.exit()          
    
    elif is_change_loc:
        if sta_code is None or lat is None or lon is None or elev is None or change_date is None:
            print("###############################################################################\n" +
                  "#### !!!Change a station location must set station-code, lat, lon, elev!!! ####\n" +
                  "###############################################################################\n")
            parser.print_help()
            sys.exit()
        else:
            inventory_gen.change_latest_location(network_code=net_code, station_code=sta_code, latitude=lat, longitude=lon, elevation=elev, created=change_date, depth=dep)
    elif is_change_ins:
        if sta_code is None or ins_key is None or change_date is None:
            print("##########################################################################\n" +
                  "#### !!!Change a station instrument must set station-code, ins_key!!! ####\n" +
                  "##########################################################################\n")
            parser.print_help()
            sys.exit()
        else:
            for instrument in InstrumentType:
                if ins_key == instrument.value[0]:
                    instrument_str = instrument.value[1]
                    inventory_gen.change_latest_instrument(network_code=net_code, station_code=sta_code, instrument_keys=instrument, changed=change_date)
                    break
            else:
                instruments_info_str2 = ""
                for instr in InstrumentType:
                    instruments_info_str2 += f'\n#### {instr.value[0]}: {instr.value[1]:<25} ####'
                print("######################################\n" +
                      "####    !!!ins-key must set!!!    ####" +
                      f"{instruments_info_str2}" +
                      "\n######################################\n")
                parser.print_help()
                sys.exit()
                
    inventory_gen.write()

    