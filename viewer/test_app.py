import subprocess
import time
import os
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Define paths
PROJECT_ROOT = '/home/jack/code/kissing_problem_viewer'
CHROMEDRIVER_PATH = '/home/jack/code/chromedriver-linux64/chromedriver'

# Start a simple HTTP server in the background
print("Starting HTTP server...")
server_process = subprocess.Popen(['python3', '-m', 'http.server', '8000'], cwd=PROJECT_ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
time.sleep(2) # Give the server a moment to start

# Configure Chrome options for headless execution
chrome_options = Options()
chrome_options.add_argument('--headless')
chrome_options.add_argument('--no-sandbox')
chrome_options.add_argument('--disable-dev-shm-usage')
chrome_options.add_argument('--window-size=1920,1080')
chrome_options.binary_location = '/home/jack/code/chrome/linux-138.0.7204.94/chrome-linux64/chrome'

# Initialize the WebDriver
service = Service(CHROMEDRIVER_PATH)
driver = None

try:
    driver = webdriver.Chrome(service=service, options=chrome_options)
    print("WebDriver initialized. Navigating to application...")

    # Navigate to the application
    driver.get('http://localhost:8000/index.html')

    # Wait for the canvas to be present
    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, "canvas"))
    )
    print("Canvas element found.")

    # Check for JavaScript errors in console logs
    print("Checking for JavaScript console errors...")
    logs = driver.get_log('browser')
    js_errors = [log for log in logs if log['level'] == 'SEVERE' and 'favicon.ico' not in log['message']]
    if js_errors:
        print("*** JavaScript Errors Found: ***")
        for error in js_errors:
            print(error)
        raise Exception("JavaScript errors detected in the browser console.")
    else:
        print("No SEVERE JavaScript errors found.")

    # Adjust zoom slider to display full website and circles
    print("Adjusting zoom slider...")
    zoom_slider = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, 'zoomSlider'))
    )
    # Set a value that should make circles visible given the large numbers in frames.json
    driver.execute_script("arguments[0].value = arguments[1]; arguments[0].dispatchEvent(new Event('input'));", zoom_slider, 0.05)
    time.sleep(1) # Give time for re-rendering

    # Take a screenshot after adjusting zoom
    screenshot_path_zoom = os.path.join(PROJECT_ROOT, 'screenshot_after_zoom.png')
    driver.save_screenshot(screenshot_path_zoom)
    print(f"Screenshot saved to {screenshot_path_zoom}")

    # Adjust frame slider and take another screenshot
    print("Adjusting frame slider...")
    frame_slider = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, 'frameSlider'))
    )
    # Move to a different frame (e.g., the last frame)
    driver.execute_script("arguments[0].value = arguments[0].max; arguments[0].dispatchEvent(new Event('input'));", frame_slider)
    time.sleep(1) # Give time for re-rendering

    screenshot_path_frame = os.path.join(PROJECT_ROOT, 'screenshot_after_frame_change.png')
    driver.save_screenshot(screenshot_path_frame)
    print(f"Screenshot saved to {screenshot_path_frame}")

    print("Selenium tests completed successfully.")

except Exception as e:
    print(f"An error occurred: {e}")
finally:
    if driver:
        print("Quitting WebDriver...")
        driver.quit()
    if server_process.poll() is None:
        print("Terminating HTTP server...")
        server_process.terminate()
    server_process.wait()
    print("HTTP server terminated.")
