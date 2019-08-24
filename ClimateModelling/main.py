"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
from user_entry import user_entry
import time


if __name__ == "__main__":
    start = time.time()
    user_entry()
    print("SUCCESS - Time to run:", time.time() - start)