# Use a base image
FROM python:3.10.13

# Set the working directory
WORKDIR /app

# Copy files
COPY . .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

RUN mkdir -p /app/outputs/imgs

# Run the application
CMD ["python", "main.py"]
