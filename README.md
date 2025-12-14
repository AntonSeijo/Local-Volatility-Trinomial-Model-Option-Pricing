# Binomial Options Pricer

Sistema completo de pricing de opciones usando el modelo binomial, implementado en C con interfaz web React.

## 🏗️ Arquitectura

- **Backend**: FastAPI + C (vía ctypes) para cálculos de alto rendimiento
- **Frontend**: React + Vite con visualización de árboles binomiales
- **Containerización**: Docker + Docker Compose

## 🚀 Inicio Rápido con Docker

### Levantar todo el stack:

```bash
docker-compose up --build
```

Esto levantará:
- **Backend** en `http://localhost:8000`
- **Frontend** en `http://localhost:80`

### Detener los servicios:

```bash
docker-compose down
```

## 💻 Desarrollo Local

### Backend

```bash
cd backend

# Compilar la librería C
gcc -shared -fPIC binomial_model.c -o binomial_model.so -lm

# Instalar dependencias Python
pip install -r requirements.txt

# Ejecutar servidor
uvicorn api:app --reload --port 8000
```

### Frontend

```bash
cd frontend

# Instalar dependencias
npm install

# Modo desarrollo
npm run dev

# Build para producción
npm run build
```

## 📡 API Endpoints

### GET `/price`
Calcula solo el precio de la opción.

**Parámetros:**
- `S0`: Precio spot
- `K`: Strike
- `T`: Tiempo hasta vencimiento (años)
- `N`: Número de pasos
- `r`: Tasa libre de riesgo
- `sigma`: Volatilidad
- `type`: `CALL` o `PUT`

**Ejemplo:**
```bash
curl "http://localhost:8000/price?S0=100&K=100&T=1&N=50&r=0.05&sigma=0.2&type=CALL"
```

### GET `/tree`
Calcula el precio y devuelve los árboles completos de precios y valores.

**Respuesta:**
```json
{
  "price": 10.4506,
  "priceTree": [[...], [...], ...],
  "valueTree": [[...], [...], ...]
}
```

## 🎨 Características del Frontend

- ✨ Diseño moderno con glassmorphism
- 📊 Visualización interactiva del árbol binomial
- 🎯 Formulario intuitivo para parámetros
- ⚡ Animaciones suaves con Framer Motion
- 📱 Diseño responsive

## 🐳 Configuración de Docker

### Variables de Entorno

**Frontend** (`.env`):
```env
VITE_API_URL=http://localhost:8000
```

Para producción, cambiar a la URL real del backend.

### Puertos

- Backend: `8000`
- Frontend: `80`

## 📝 Notas

- El visualizador de árbol está optimizado para N ≤ 15
- Para N > 15, se recomienda usar solo el endpoint `/price`
- El backend usa C compilado para máximo rendimiento

## 🛠️ Tecnologías

- **Backend**: Python 3.11, FastAPI, ctypes, NumPy
- **Frontend**: React 18, Vite, Axios, Recharts, Framer Motion
- **Cálculo**: C (GCC)
- **Deployment**: Docker, Nginx
