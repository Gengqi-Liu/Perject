import nest_asyncio
nest_asyncio.apply()

from bleak import BleakClient
import asyncio
import keyboard
import time
import math
import logging
#import platform
from bleak import _logger as logger
import socket

# Setup parameters
CHARACTERISTIC_UUID = ("6e400003-b5a3-f393-e0a9-e50e24dcca9e")
ADDRESS             = ("F9:5D:57:CE:6F:68") # Bluefruit 52
# ADDRESS             = ("F8:AA:DF:51:6F:F2") # Headtracker_3 Nachbau
# ADDRESS             = ("E0:E5:D2:19:3C:61") # Bluefruit 52_2, Papa
# ADDRESS             = ("FE:C4:FE:10:50:0E") # Joels Headtracker
# ADDRESS             = ("F2:64:83:4A:0A:DA") # HT_BNO085


uuid_battery_service = '0000180f-0000-1000-8000-00805f9b34fb'
uuid_battery_level_characteristic = '00002a19-0000-1000-8000-00805f9b34fb'

fAngleRefHor        = 0.0
fAngleRefVer        = 0.0
fDuration           = 60*120
PI                  = 3.141592653589793
time1               = 0

#global sock, UDP_IP, UDP_PORT
UDP_IP = "127.0.0.1"
UDP_PORT = 5005
print("UDP target IP: %s" % UDP_IP)
print("UDP target port: %s" % UDP_PORT)
sock = socket.socket(socket.AF_INET, # Internet
                     socket.SOCK_DGRAM) # UDP

# Notification function
#def notification_handler(self, sender: str, data:any):
def notification_handler(sender,data):

    global fAngleRefHor, fAngleRefVer
    global time1,time2
    global sock, UDP_IP, UDP_PORT
    time2       = time.perf_counter()*1000
    timediff    = time2-time1
    #print(timediff)
    time1       = time2
    if data!=b'\x00\x00\x00\x00c':
        data = int.from_bytes(data, byteorder='little', signed=False)
        fAngleVer0  = math.floor(data/pow(2,13))
        fAngleHor0  = data-fAngleVer0*pow(2,13)   
        fAngleHor0  = 180/PI*fAngleHor0/1000
        fAngleVer0  = 180/PI*fAngleVer0/100
        fAngleHor   = (fAngleHor0-fAngleRefHor-180) % 360 - 180;
        fAngleVer   = (fAngleVer0-fAngleRefVer-180) % 360 - 180;   
        fAngleHor   = '%+8.3f' % fAngleHor
        fAngleVer   = '%+8.3f' % fAngleVer
        bAngleHor   = bytes(fAngleHor, 'utf-8')
        bAngleVer   = bytes(fAngleVer, 'utf-8')
        sock.sendto(bAngleHor+bAngleVer, (UDP_IP, UDP_PORT))
    if keyboard.is_pressed('alt+c'):
        fAngleRefHor = fAngleHor0
        fAngleRefVer = fAngleVer0
    if keyboard.is_pressed('alt+1'):
        f = open("FilterNumber.txt", "w")
        f.write("1")
        f.close()
    if keyboard.is_pressed('alt+2'):
        f = open("FilterNumber.txt", "w")
        f.write("2")
        f.close()

# Loop
async def run(address,debug=False):
    if debug:
        import sys
        l = logging.getLogger("asyncio")
        l.setLevel(logging.DEBUG)
        h = logging.StreamHandler(sys.stdout)
        h.setLevel(logging.DEBUG)
        l.addHandler(h)
        logger.addHandler(h)
    async with BleakClient(address) as client:
        global fAngleRefHor, fAngleRefVer
        await client.is_connected()
        print('Connected to headtracker ...')

        print('Try to read battery level')
        try:
            battery_level = await client.read_gatt_char(uuid_battery_level_characteristic)
            battery_level = int.from_bytes(battery_level, byteorder='big')
            print('Battery level: ',str(battery_level),'%')
            print('Battery level read')
        except:
            print('Battery level could not be read')

        bInitAngles = False
        while not bInitAngles:
            print('Try to read initial values ...')
            data = await client.read_gatt_char(CHARACTERISTIC_UUID)
            if data!=b'\x00\x00\x00\x00':
                data = int.from_bytes(data, byteorder='little', signed=False)
                fAngleRefVer    = math.floor(data/pow(2,13))
                fAngleRefHor    = data-fAngleRefVer*pow(2,13)
                fAngleRefHor    = 180/PI*fAngleRefHor/1000
                fAngleRefVer    = 180/PI*fAngleRefVer/100
                bInitAngles     = True
                print('Initial values read ...')
        await client.start_notify(CHARACTERISTIC_UUID, notification_handler)
        print('Headtracker notification started ...')
        print('Press alt+c to center ...')
        await asyncio.sleep(fDuration, loop=loop)
        await client.stop_notify(CHARACTERISTIC_UUID)
        print('Headtracker notification stopped ...')
        print("End headtracker session!")
        
       
# Main function
if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    loop.set_debug(True)
    loop.run_until_complete(run(ADDRESS,True))
    loop.close()