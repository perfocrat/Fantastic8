const geometryDefinitions = {
    rectangle: {
        name: "Rectangle",
        dimensionUnit: "LENGTH · nm",
        dimensions: [
            { key: "width", label: "Width", value: 8.0, unit: "nm" },
            { key: "height", label: "Height", value: 5.0, unit: "nm" }
        ]
    },
    triangle: {
        name: "Equilateral Triangle",
        dimensionUnit: "LENGTH · nm",
        dimensions: [{ key: "side", label: "Side length", value: 7.0, unit: "nm" }]
    },
    circle: {
        name: "Circle",
        dimensionUnit: "LENGTH · nm",
        dimensions: [{ key: "radius", label: "Radius", value: 4.0, unit: "nm" }]
    },
    sector: {
        name: "Circular Sector",
        dimensionUnit: "LENGTH / ANGLE",
        dimensions: [
            { key: "radius", label: "Radius", value: 4.0, unit: "nm" },
            { key: "angle", label: "Sector angle", value: 90, unit: "°" }
        ]
    },
    annulus: {
        name: "Annulus",
        dimensionUnit: "LENGTH · nm",
        dimensions: [
            { key: "innerRadius", label: "Inner radius", value: 2.0, unit: "nm" },
            { key: "outerRadius", label: "Outer radius", value: 5.0, unit: "nm" }
        ]
    },
    ellipse: {
        name: "Ellipse",
        dimensionUnit: "LENGTH · nm",
        dimensions: [
            { key: "semiMajor", label: "Semi-major axis", value: 5.0, unit: "nm" },
            { key: "semiMinor", label: "Semi-minor axis", value: 3.0, unit: "nm" }
        ]
    }
};

const demoResults = {
    rectangle: { energy: "12.4", nextEnergy: "18.7" },
    triangle: { energy: "14.2", nextEnergy: "21.6" },
    circle: { energy: "10.8", nextEnergy: "16.3" },
    sector: { energy: "11.6", nextEnergy: "17.9" },
    annulus: { energy: "13.1", nextEnergy: "19.4" },
    ellipse: { energy: "9.7", nextEnergy: "15.2" }
};

const state = {
    geometry: "rectangle",
    quantumNumbers: { n: 1, m: 2 },
    dimensions: Object.fromEntries(
        Object.entries(geometryDefinitions).map(([key, definition]) => [
            key,
            Object.fromEntries(definition.dimensions.map(({ key: dimensionKey, value }) => [dimensionKey, value]))
        ])
    )
};

const elements = {
    geometryCards: document.querySelectorAll(".geometry-card"),
    geometryOptions: document.querySelectorAll(".geometry-option"),
    brandHome: document.getElementById("brand-home"),
    landing: document.getElementById("geometry-landing"),
    workspace: document.getElementById("workspace"),
    backButton: document.getElementById("back-to-geometries"),
    workspaceTitle: document.getElementById("workspace-title"),
    geometryName: document.getElementById("active-geometry-name"),
    dimensionUnit: document.getElementById("dimension-unit"),
    dimensionControls: document.getElementById("dimension-controls"),
    nSlider: document.getElementById("quantum-n"),
    mSlider: document.getElementById("quantum-m"),
    nValue: document.getElementById("n-value"),
    mValue: document.getElementById("m-value"),
    nMaximumLabel: document.getElementById("n-maximum-label"),
    mMaximumLabel: document.getElementById("m-maximum-label"),
    nLabelContext: document.getElementById("n-label-context"),
    mLabelContext: document.getElementById("m-label-context"),
    constraintNote: document.getElementById("constraint-note"),
    energyValue: document.getElementById("energy-value"),
    nextEnergyValue: document.getElementById("next-energy-value")
};

function isTriangleSelected() {
    return state.geometry === "triangle";
}

function renderDimensions() {
    const definition = geometryDefinitions[state.geometry];
    elements.dimensionUnit.textContent = definition.dimensionUnit;
    elements.dimensionControls.replaceChildren();

    definition.dimensions.forEach(({ key, label, unit }) => {
        const control = document.createElement("div");
        control.className = "dimension-control";

        const inputId = `dimension-${key}`;
        const fieldLabel = document.createElement("label");
        fieldLabel.htmlFor = inputId;
        fieldLabel.textContent = label;

        const inputWrap = document.createElement("div");
        inputWrap.className = "dimension-input-wrap";

        const input = document.createElement("input");
        input.className = "dimension-input";
        input.id = inputId;
        input.name = key;
        input.type = "number";
        input.min = "0";
        input.step = unit === "°" ? "1" : "0.1";
        input.value = state.dimensions[state.geometry][key];
        input.setAttribute("inputmode", "decimal");
        input.addEventListener("input", () => {
            if (input.value !== "") {
                state.dimensions[state.geometry][key] = Number(input.value);
            }
        });
        input.addEventListener("change", refreshPreview);

        const unitLabel = document.createElement("span");
        unitLabel.className = "dimension-unit";
        unitLabel.textContent = unit;

        inputWrap.append(input, unitLabel);
        control.append(fieldLabel, inputWrap);
        elements.dimensionControls.append(control);
    });
}

function renderQuantumNumbers() {
    const triangleSelected = isTriangleSelected();
    const nMaximum = triangleSelected ? state.quantumNumbers.m : 10;
    const mMaximum = triangleSelected ? 5 : 10;

    elements.mSlider.min = "1";
    elements.mSlider.max = String(mMaximum);
    elements.mSlider.value = String(state.quantumNumbers.m);
    elements.nSlider.min = "1";
    elements.nSlider.max = String(nMaximum);
    elements.nSlider.value = String(state.quantumNumbers.n);
    elements.nValue.value = String(state.quantumNumbers.n);
    elements.nValue.textContent = String(state.quantumNumbers.n);
    elements.mValue.value = String(state.quantumNumbers.m);
    elements.mValue.textContent = String(state.quantumNumbers.m);
    elements.nMaximumLabel.textContent = String(nMaximum).padStart(2, "0");
    elements.mMaximumLabel.textContent = String(mMaximum).padStart(2, "0");
    elements.nLabelContext.textContent = triangleSelected ? "triangle index" : "principal index";
    elements.mLabelContext.textContent = triangleSelected ? "upper index" : "secondary index";
    elements.constraintNote.hidden = !triangleSelected;
}

function renderEnergyResults(result) {
    elements.energyValue.innerHTML = `${result.values.energy} <small>meV</small>`;
    elements.nextEnergyValue.innerHTML = `${result.values.nextEnergy} <small>meV</small>`;
}

function renderGeometry() {
    const definition = geometryDefinitions[state.geometry];
    elements.geometryOptions.forEach((option) => {
        const selected = option.dataset.geometry === state.geometry;
        option.classList.toggle("is-active", selected);
        option.setAttribute("aria-pressed", String(selected));
    });

    elements.workspaceTitle.textContent = definition.name;
    elements.geometryName.textContent = definition.name;
    renderDimensions();
    renderQuantumNumbers();
}

function selectGeometry(geometry) {
    if (!geometryDefinitions[geometry]) return;

    state.geometry = geometry;
    if (isTriangleSelected()) {
        state.quantumNumbers.m = Math.min(state.quantumNumbers.m, 5);
        state.quantumNumbers.n = Math.min(state.quantumNumbers.n, state.quantumNumbers.m);
    }

    renderGeometry();
    showWorkspace();
    refreshPreview();
}

function showWorkspace() {
    elements.landing.hidden = true;
    elements.workspace.hidden = false;
    elements.workspaceTitle.focus?.({ preventScroll: true });
    window.scrollTo({ top: 0, behavior: "smooth" });
}

function showLanding() {
    elements.workspace.hidden = true;
    elements.landing.hidden = false;
    window.scrollTo({ top: 0, behavior: "smooth" });
    elements.landing.querySelector(".geometry-card")?.focus({ preventScroll: true });
}

function onMChange(event) {
    state.quantumNumbers.m = Number(event.currentTarget.value);
    if (isTriangleSelected() && state.quantumNumbers.n > state.quantumNumbers.m) {
        state.quantumNumbers.n = state.quantumNumbers.m;
    }
    renderQuantumNumbers();
    refreshPreview();
}

function onNChange(event) {
    state.quantumNumbers.n = Number(event.currentTarget.value);
    renderQuantumNumbers();
    refreshPreview();
}

function buildRequestPayload() {
    return {
        geometry: state.geometry,
        dimensions: { ...state.dimensions[state.geometry] },
        quantumNumbers: { ...state.quantumNumbers }
    };
}

async function requestPreview(payload) {
    // Replace this demo provider with fetch() when a backend endpoint is available.
    return Promise.resolve({
        source: "demo",
        geometry: payload.geometry,
        values: demoResults[payload.geometry]
    });
}

async function refreshPreview() {
    const requestPayload = buildRequestPayload();
    const result = await requestPreview(requestPayload);
    renderEnergyResults(result);
}

elements.geometryOptions.forEach((option) => {
    option.addEventListener("click", () => selectGeometry(option.dataset.geometry));
});
elements.geometryCards.forEach((card) => {
    card.addEventListener("click", () => selectGeometry(card.dataset.geometry));
});
elements.backButton.addEventListener("click", showLanding);
elements.brandHome.addEventListener("click", (event) => {
    if (!elements.workspace.hidden) {
        event.preventDefault();
        showLanding();
    }
});
elements.mSlider.addEventListener("input", onMChange);
elements.nSlider.addEventListener("input", onNChange);
document.getElementById("parameter-form").addEventListener("submit", (event) => event.preventDefault());

renderGeometry();