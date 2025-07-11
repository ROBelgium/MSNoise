"""
MSNoise Web Administration Interface - Consolidated Version

This module provides a Flask-Admin based web interface for managing MSNoise
configurations, workflow chains, filters, and data processing jobs.

Key Features:
- Unified Configuration Management (global and workflow-specific)
- Dynamic Workflow Chain Creation and Management
- Enhanced Filter Management with Workflow Associations
- Simplified Base Classes with Modern Flask-Admin
- Generic Workflow Handling driven by Database Schema
"""

import datetime
import logging
import os
import sys
import traceback
from collections import OrderedDict

from flask import Flask, request, redirect, url_for, flash, jsonify, render_template_string
from flask_admin import Admin, BaseView, expose, AdminIndexView
from flask_admin.contrib.sqla import ModelView
from flask_admin.form import Select2Widget
from flask_admin.model import typefmt
from flask_admin import helpers as admin_helpers
from flask_babel import Babel, gettext, ngettext, lazy_gettext
from markupsafe import Markup
SafeMarkup = Markup
from sqlalchemy import and_, or_, func, desc, asc
from sqlalchemy.orm import sessionmaker
from wtforms import Form, StringField, TextAreaField, SelectField, BooleanField, IntegerField, FloatField
from wtforms.validators import DataRequired, Optional, NumberRange, ValidationError
from wtforms.widgets import TextArea, Select, CheckboxInput, NumberInput

from .api import connect, get_config, get_logger
from .api import get_config_categories_definition
from .msnoise_table_def import declare_tables


# Initialize Flask app and extensions
app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'msnoise-admin-secret-key')
app.config['FLASK_ADMIN_SWATCH'] = 'cerulean'

# Babel for internationalization
babel = Babel(app)

# Database connection
db_session = connect()

# Get schema from declare_tables
schema = declare_tables()

# Logger
logger = get_logger('msnoise.admin', loglevel="INFO")


# ============================================================================
# Base Classes and Utilities
# ============================================================================

class BaseModelView(ModelView):
    """
    Enhanced base model view with common functionality for all MSNoise models.
    Simplified from the original complex implementation.
    """

    # Common settings
    can_export = True
    can_view_details = True
    page_size = 50
    can_set_page_size = True

    # Default column formatters
    column_formatters = {
        'used': lambda view, context, model, name: SafeMarkup('✓' if model.used else '✗') if hasattr(model, 'used') else '',
    }

    def __init__(self, model, session, *args, **kwargs):
        super(BaseModelView, self).__init__(model, session, *args, **kwargs)
        self.session = session

    def get_query(self):
        """Override to add common filters"""
        return self.session.query(self.model)

    def get_count_query(self):
        """Override to add common filters"""
        return self.session.query(func.count('*')).select_from(self.model)


class MSNoiseAdminIndexView(AdminIndexView):
    """Custom admin index view with MSNoise-specific dashboard"""

    @expose('/')
    def index(self):
        """Display MSNoise dashboard"""
        # Get basic statistics
        stats = {}
        try:
            stats['stations'] = db_session.query(schema.Station).filter(schema.Station.used == True).count()
            # stats['filters'] = db_session.query(schema.Filter).filter(schema.Filter.used == True).count()
            # stats['workflow_steps'] = db_session.query(schema.WorkflowSteps).filter(schema.WorkflowSteps.used == True).count()
            stats['jobs_total'] = db_session.query(schema.Job).count()
            stats['jobs_todo'] = db_session.query(schema.Job).filter(schema.Job.flag == 'T').count()
            stats['jobs_done'] = db_session.query(schema.Job).filter(schema.Job.flag == 'D').count()
            stats['config_params'] = db_session.query(schema.Config).count()
        except Exception as e:
            logger.error(f"Error getting dashboard stats: {e}")
            stats = {}

        return self.render('admin/index.html', stats=stats)


# ============================================================================
# 1. Single ConfigView - Unified Configuration Management
# ============================================================================

class ConfigView(BaseModelView):
    """
    Unified configuration management for both global and workflow-specific settings.
    This replaces all the separate DvvMwcsView, DvvStretchingView, etc. classes.
    """

    column_list = ['category', 'name', 'set_number', 'value', 'description']
    column_searchable_list = ['name', 'category', 'value', 'description']
    column_filters = ['category', 'set_number']
    column_sortable_list = ['name', 'category', 'set_number']

    form_columns = ['name', 'category', 'set_number', 'value', 'param_type',
                    'default_value', 'description', 'units', 'possible_values']

    column_descriptions = {
        'name': 'Parameter name (e.g., maxlag, dtt_minlag)',
        'category': 'Configuration category (global, mwcs, stretching, etc.)',
        'set_number': 'Configuration set number (NULL for global configs)',
        'value': 'Current parameter value',
        'param_type': 'Data type of the parameter',
        'description': 'Description of what this parameter does',
    }

    form_choices = {
        'category': get_config_categories_definition(),
        'param_type': [
            ('str', 'String'),
            ('int', 'Integer'),
            ('float', 'Float'),
            ('bool', 'Boolean'),
        ]
    }

    # Custom column formatters
    column_formatters = dict(BaseModelView.column_formatters, **{
        'value': lambda view, context, model, name: SafeMarkup(f'<code>{model.value}</code>') if model.value else '',
        'category': lambda view, context, model, name: SafeMarkup(f'<span class="label label-info">{model.category}</span>'),
        'param_type': lambda view, context, model, name: SafeMarkup(f'<span class="label label-default">{model.param_type}</span>'),
        'set_number': lambda view, context, model, name: str(model.set_number) if model.set_number is not None else 'Global',
    })

    def __init__(self, session, *args, **kwargs):
        super(ConfigView, self).__init__(schema.Config, session, *args, **kwargs)

    def create_model(self, form):
        """Create new configuration parameter with validation"""
        try:
            model = self.model()
            form.populate_obj(model)

            # Validate value type
            if model.param_type == 'int':
                int(model.value)
            elif model.param_type == 'float':
                float(model.value)
            elif model.param_type == 'bool':
                model.value = str(model.value).lower() in ('true', 'yes', '1', 'on')

            # Validate against possible values if defined
            model.validate_value()

            self.session.add(model)
            self.session.commit()
            return model
        except Exception as e:
            self.session.rollback()
            flash(f'Failed to create configuration: {str(e)}', 'error')
            return False

    def update_model(self, form, model):
        """Update configuration parameter with validation"""
        try:
            form.populate_obj(model)

            # Validate value type
            if model.param_type == 'int':
                int(model.value)
            elif model.param_type == 'float':
                float(model.value)
            elif model.param_type == 'bool':
                model.value = str(model.value).lower() in ('true', 'yes', '1', 'on')

            # Validate against possible values if defined
            model.validate_value()

            self.session.commit()
            return True
        except Exception as e:
            self.session.rollback()
            flash(f'Failed to update configuration: {str(e)}', 'error')
            return False

    @expose('/bulk_update/', methods=['GET', 'POST'])
    def bulk_update(self):
        """Bulk update configuration parameters"""
        if request.method == 'POST':
            try:
                updates = request.form.to_dict()
                for key, value in updates.items():
                    if key.startswith('config_'):
                        config_id = int(key.split('_')[1])
                        config = self.session.query(schema.Config).get(config_id)
                        if config:
                            config.value = value
                            config.updated_at = datetime.datetime.utcnow()

                self.session.commit()
                flash('Configuration updated successfully', 'success')
            except Exception as e:
                self.session.rollback()
                flash(f'Failed to update configuration: {str(e)}', 'error')

            return redirect(url_for('.index_view'))

        # GET request - show bulk update form
        configs = self.session.query(schema.Config).all()
        return self.render('admin/config_bulk_update.html', configs=configs)


# ============================================================================
# 4. Simplified Base Classes - Standard Model Views
# ============================================================================

class StationView(BaseModelView):
    """Simplified station management"""

    column_list = ['ref', 'net', 'sta', 'X', 'Y', 'altitude', 'coordinates', 'used']
    column_searchable_list = ['net', 'sta']
    column_filters = ['net', 'coordinates', 'used']
    column_sortable_list = ['ref', 'net', 'sta', 'X', 'Y', 'used']

    form_columns = ['net', 'sta', 'X', 'Y', 'altitude', 'coordinates',
                    'used_location_codes', 'used_channel_names', 'used']

    def __init__(self, session, *args, **kwargs):
        super(StationView, self).__init__(schema.Station, session, *args, **kwargs)


class JobView(BaseModelView):
    """Simplified job management"""

    column_list = ['ref', 'day', 'pair', 'jobtype', 'flag', 'lastmod']
    column_searchable_list = ['day', 'pair', 'jobtype']
    column_filters = ['jobtype', 'flag']
    column_sortable_list = ['ref', 'day', 'pair', 'jobtype', 'flag', 'lastmod']

    form_columns = ['day', 'pair', 'jobtype', 'flag']

    # Custom column formatters
    column_formatters = dict(BaseModelView.column_formatters, **{
        'flag': lambda view, context, model, name: SafeMarkup(
            f'<span class="label label-{"success" if model.flag == "D" else "warning" if model.flag == "I" else "default"}">'
            f'{"Done" if model.flag == "D" else "In Progress" if model.flag == "I" else "Todo"}</span>'
        ),
        'lastmod': lambda view, context, model, name: model.lastmod.strftime('%Y-%m-%d %H:%M:%S') if model.lastmod else '',
    })

    def __init__(self, session, *args, **kwargs):
        super(JobView, self).__init__(schema.Job, session, *args, **kwargs)


class DataAvailabilityView(BaseModelView):
    """Simplified data availability management"""

    column_list = ['ref', 'net', 'sta', 'loc', 'chan', 'starttime', 'endtime', 'samplerate', 'flag']
    column_searchable_list = ['net', 'sta', 'loc', 'chan', 'file']
    column_filters = ['net', 'sta', 'chan', 'flag']
    column_sortable_list = ['ref', 'net', 'sta', 'starttime', 'endtime', 'samplerate']

    # Read-only view
    can_create = False
    can_edit = False
    can_delete = False

    def __init__(self, session, *args, **kwargs):
        super(DataAvailabilityView, self).__init__(schema.DataAvailability, session, *args, **kwargs)

from flask_admin import BaseView, expose
from flask_admin.form import BaseForm
from flask import request, redirect, url_for, flash
from wtforms import Form, StringField, FloatField, IntegerField, BooleanField, SelectField, TextAreaField
from wtforms.validators import DataRequired, Optional
import datetime

class ConfigSetView(BaseView):
    """
    View for managing configuration sets as single forms.
    Allows editing all parameters in a set together.
    """

    def __init__(self, session, *args, **kwargs):
        super(ConfigSetView, self).__init__(*args, **kwargs)
        self.session = session

    @expose('/')
    def index(self):
        """Display all configuration sets grouped by category"""

        """Display all configuration sets grouped by category"""
        from .api import get_config_sets_organized

        ordered_categories = get_config_sets_organized(self.session)

        return self.render('admin/config_sets_index.html',
                           categories=ordered_categories)


    @expose('/create/<category>/', methods=['GET', 'POST'])
    def create_set(self, category):
        """Create a new configuration set for a category"""

        from .api import create_config_set
        set_number = create_config_set(self.session, category)


        if set_number is not None:
            self.session.commit()
            flash(f'Configuration set "{category}" #{set_number} created successfully', 'success')
            return redirect(url_for('.edit_set', category=category, set_number=set_number))

        return

    @expose('/view/<category>/<int:set_number>/')
    def view_set(self, category, set_number):
        """View configuration set in read-only mode"""

        # Get all configs for this category and set
        configs = self.session.query(schema.Config).filter(
            schema.Config.category == category,
            schema.Config.set_number == set_number
        ).order_by(schema.Config.name).all()

        if not configs:
            flash(f'Configuration set "{category}" #{set_number} not found', 'error')
            return redirect(url_for('.index'))

        # Group configs by logical sections if needed
        config_groups = {}
        for config in configs:
            config.description = Markup(config.description)
            # You can add logic here to group configs by prefix or type
            group_name = "Configuration Parameters"
            if group_name not in config_groups:
                config_groups[group_name] = []

            config_groups[group_name].append(config)

        return self.render('admin/config_set_view.html',
                           category=category,
                           set_number=set_number,
                           configs=configs,
                           config_groups=config_groups)


    @expose('/edit/<category>/<int:set_number>/', methods=['GET', 'POST'])
    def edit_set(self, category, set_number):
        """Edit a complete configuration set as a single form"""

        # Get all config parameters for this set
        configs = self.session.query(schema.Config).filter(
            schema.Config.category == category,
            schema.Config.set_number == set_number
        ).order_by(schema.Config.name).all()

        if not configs:
            flash(f'Configuration set "{category}" #{set_number} not found', 'error')
            return redirect(url_for('.index'))

        # Create dynamic form class
        form_class = self._create_dynamic_form(configs)
        # Initialize form with request data for POST requests
        if request.method == 'POST':
            form = form_class(request.form)
        else:
            form = form_class()

        if request.method == 'POST' and form.validate():
            try:
                # Update all config values
                for config in configs:
                    field_name = f"param_{config.ref}"
                    if hasattr(form, field_name):
                        field = getattr(form, field_name)
                        new_value = field.data

                        # Handle different field types and preserve existing values if new value is empty
                        if config.param_type == 'bool':
                            # For boolean fields
                            if isinstance(field, BooleanField):
                                # WTForms BooleanField returns True/False
                                new_value = 'Y' if new_value else 'N'
                            elif isinstance(field, SelectField):
                                # SelectField for boolean with Y/N choices
                                new_value = str(new_value) if new_value else config.value
                            else:
                                new_value = str(new_value) if new_value else config.value
                        elif config.param_type == 'int':
                            if new_value is not None and str(new_value).strip():
                                try:
                                    new_value = str(int(new_value))
                                except (ValueError, TypeError):
                                    # Keep original value if conversion fails
                                    new_value = config.value
                            else:
                                # Keep original value if empty
                                new_value = config.value
                        elif config.param_type == 'float':
                            if new_value is not None and str(new_value).strip():
                                try:
                                    new_value = str(float(new_value))
                                except (ValueError, TypeError):
                                    # Keep original value if conversion fails
                                    new_value = config.value
                            else:
                                # Keep original value if empty
                                new_value = config.value
                        else:  # string type
                            if new_value is not None:
                                new_value = str(new_value).strip()
                                # Allow empty strings for string fields
                            else:
                                new_value = config.value

                        # Only update if we have a valid new value
                        if new_value is not None:
                            config.value = new_value
                            config.updated_at = datetime.datetime.utcnow()

                            # Validate against possible values (only for non-empty values)
                            if new_value.strip():
                                try:
                                    config.validate_value()
                                except ValueError as ve:
                                    flash(f'Validation error for {config.name}: {str(ve)}', 'error')
                                    raise ve
                            # flash(f'Updated parameter "{config.name}"="{new_value}"', 'success')
                            self.session.commit()
                self.session.commit()
                flash(f'Configuration set "{category}" #{set_number} updated successfully', 'success')
                return redirect(url_for('.view_set', category=category, set_number=set_number))

            except Exception as e:
                self.session.rollback()
                flash(f'Error updating configuration: {str(e)}', 'error')

        else:
            # Populate form with current values (GET request or validation failed)
            for config in configs:
                config.description = Markup(config.description)
                field_name = f"param_{config.ref}"
                if hasattr(form, field_name):
                    field = getattr(form, field_name)
                    current_value = config.value if config.value is not None else ''
                    field.description = config.description
                    if config.param_type == 'bool':
                        # Handle boolean conversion for different formats
                        if isinstance(field, BooleanField):
                            field.data = current_value.upper() in ('Y', 'TRUE', 'YES', '1', 'ON')
                        else:  # SelectField for boolean
                            field.data = current_value
                    elif config.param_type in ['int', 'float']:
                        try:
                            if current_value and str(current_value).strip():
                                if config.param_type == 'float':
                                    field.data = float(current_value)
                                else:
                                    field.data = int(current_value)
                            else:
                                field.data = None
                        except (ValueError, TypeError):
                            # If conversion fails, set as string for display
                            field.data = current_value
                    else:
                        field.data = current_value

            # Show validation errors if any
            if request.method == 'POST':
                for field_name, errors in form.errors.items():
                    for error in errors:
                        flash(f'Validation error in {field_name}: {error}', 'error')

        return self.render('admin/config_set_edit.html',
                           form=form,
                           configs=configs,
                           category=category,
                           set_number=set_number)

    @expose('/copy/<category>/<int:set_number>/', methods=['GET', 'POST'])
    def copy_set(self, category, set_number):
        """Copy a configuration set to create a new one"""

        if request.method == 'POST':
            new_category = request.form.get('new_category', category)
            new_set_number = int(request.form.get('new_set_number', set_number + 1))

            try:
                # Check if target already exists
                existing = self.session.query(schema.Config).filter(
                    schema.Config.category == new_category,
                    schema.Config.set_number == new_set_number
                ).first()

                if existing:
                    flash(f'Configuration set "{new_category}" #{new_set_number} already exists', 'error')
                    return redirect(url_for('.index'))

                # Get source configs
                source_configs = self.session.query(schema.Config).filter(
                    schema.Config.category == category,
                    schema.Config.set_number == set_number
                ).all()

                if not source_configs:
                    flash(f'Source configuration set not found', 'error')
                    return redirect(url_for('.index'))

                # Copy configs
                for config in source_configs:
                    new_config = schema.Config(
                        name=config.name,
                        category=new_category,
                        set_number=new_set_number,
                        value=config.value,
                        param_type=config.param_type,
                        default_value=config.default_value,
                        description=config.description,
                        units=config.units,
                        possible_values=config.possible_values,
                    )
                    self.session.add(new_config)

                self.session.commit()
                flash(f'Configuration set copied to "{new_category}" #{new_set_number}', 'success')
                return redirect(url_for('.edit_set', category=new_category, set_number=new_set_number))

            except Exception as e:
                self.session.rollback()
                flash(f'Error copying configuration: {str(e)}', 'error')

        return self.render('admin/config_set_copy.html',
                           category=category,
                           set_number=set_number)

    @expose('/delete/<category>/<int:set_number>/', methods=['POST'])
    def delete_set(self, category, set_number):
        """Delete a configuration set"""

        try:
            deleted_count = self.session.query(schema.Config).filter(
                schema.Config.category == category,
                schema.Config.set_number == set_number
            ).delete()
            self.session.commit()
            if deleted_count > 0:
                self.session.commit()
                flash(f'Configuration set "{category}" #{set_number} deleted ({deleted_count} parameters)', 'success')
            else:
                flash(f'Configuration set "{category}" #{set_number} not found', 'error')

        except Exception as e:
            self.session.rollback()
            flash(f'Error deleting configuration: {str(e)}', 'error')

        return redirect(url_for('.index'))

    def _create_dynamic_form(self, configs):
        """Create a dynamic WTForms class based on config parameters"""

        class DynamicConfigForm(Form):
            pass

        for config in configs:
            field_name = f"param_{config.ref}"

            # Determine field type based on param_type
            if config.param_type == 'int':
                field = IntegerField(
                    label=config.name,
                    description=config.description,
                    validators=[Optional()]  # Make all fields optional to prevent erasure
                )
            elif config.param_type == 'float':
                field = FloatField(
                    label=config.name,
                    description=config.description,
                    validators=[Optional()]  # Make all fields optional to prevent erasure
                )
            elif config.param_type == 'bool':
                # Check if we should use a dropdown or checkbox
                field = BooleanField(
                    label=config.name,
                    description=config.description
                )
            else:  # string type
                # Check if there are possible values (dropdown)
                if config.possible_values:
                    # Parse possible values and ensure current value is included
                    possible_vals = [v.strip() for v in config.possible_values.split('/') if v.strip()]

                    # Make sure the current value is in the choices (for validation)
                    current_value = config.value.strip() if config.value else ''
                    if current_value and current_value not in possible_vals:
                        possible_vals.append(current_value)

                    # Create choices as (value, label) tuples
                    choices = [(v, v) for v in possible_vals]

                    # Add empty choice for optional fields
                    if not current_value or 'optional' in (config.description or '').lower():
                        choices.insert(0, ('', '-- Select --'))

                    field = SelectField(
                        label=config.name,
                        description=config.description,
                        choices=choices,
                        validators=[Optional()]
                    )
                else:
                    # Long text or regular string
                    if config.description and len(config.description) > 100:
                        field = TextAreaField(
                            label=config.name,
                            description=config.description,
                            validators=[Optional()]
                        )
                    else:
                        field = StringField(
                            label=config.name,
                            description=config.description,
                            validators=[Optional()]
                        )

            # # Add units to description if available
            # if config.units:
            #     current_desc = field.description or ""
            #     field.description = f"{current_desc} (Units: {config.units})".strip()
            #
            # # Add possible values to description if available and not already a select field
            # if config.possible_values and not isinstance(field, SelectField):
            #     current_desc = field.description or ""
            #     field.description = f"{current_desc} Possible values: {config.possible_values}".strip()

            setattr(DynamicConfigForm, field_name, field)

        return DynamicConfigForm



CONFIG_BULK_UPDATE_TEMPLATE = """
{% extends 'admin/base.html' %}

{% block body %}
<div class="container-fluid">
    <h2>Bulk Configuration Update</h2>
    
    <form method="post">
        {% for config in configs %}
        <div class="form-group">
            <label>{{ config.name }} ({{ config.category }})</label>
            <input type="text" name="config_{{ config.ref }}" value="{{ config.value }}" class="form-control">
            <small class="help-block">{{ config.description }}</small>
        </div>
        {% endfor %}
        
        <button type="submit" class="btn btn-primary">Update All</button>
    </form>
</div>
{% endblock %}
"""

MSNOISE_INDEX_TEMPLATE = """
{% extends 'admin/index.html' %}

{% block body %}
<div class="container-fluid">
    <h2>MSNoise Dashboard</h2>
    
    <div class="row">
        <div class="col-md-3">
            <div class="panel panel-primary">
                <div class="panel-heading">Stations</div>
                <div class="panel-body">
                    <h3>{{ stats.get('stations', 0) }}</h3>
                    <p>Active stations</p>
                </div>
            </div>
        </div>
       
        
        <div class="col-md-3">
            <div class="panel panel-success">
                <div class="panel-heading">Workflows</div>
                <div class="panel-body">
                    <h3>{{ stats.get('workflow_steps', 0) }}</h3>
                    <p>Active workflow steps</p>
                </div>
            </div>
        </div>
        
        <div class="col-md-3">
            <div class="panel panel-warning">
                <div class="panel-heading">Jobs</div>
                <div class="panel-body">
                    <h3>{{ stats.get('jobs_todo', 0) }} / {{ stats.get('jobs_total', 0) }}</h3>
                    <p>Todo / Total jobs</p>
                </div>
            </div>
        </div>
    </div>
</div>
{% endblock %}
"""

import csv
import os
from flask import request, redirect, url_for, flash, render_template_string, render_template
from flask_admin import BaseView, expose
from wtforms import Form, SelectField, StringField, TextAreaField, IntegerField, validators

# ============================================================================
# Main Application Setup
# ============================================================================

def create_admin_app():
    """Create and configure the Flask-Admin application"""

    # Initialize Flask-Admin
    admin = Admin(
        app,
        name='MSNoise Admin',
        template_mode='bootstrap3',
        index_view=MSNoiseAdminIndexView(),
        base_template='admin/base.html'
    )

    # Add model views
    # admin.add_view(FilterView(db_session, name='Filters', category='Configuration'))

    admin.add_view(ConfigSetView(db_session, name='Configuration Sets', endpoint='config_sets', category='Configuration'))
    admin.add_view(StationView(db_session, name='Stations', category='Configuration'))
    admin.add_view(DataAvailabilityView(db_session, name='Data Availability', category='Data'))

    admin.add_view(JobView(db_session, name='Jobs', category='Processing'))

    admin.add_view(ConfigView(db_session, name='Configuration (Expert)', category='Expert'))
    # Add API view
    # admin.add_view(WorkflowAPIView(name='Workflow API', category='API'))

    # Add custom templates
    app.jinja_env.globals['CONFIG_BULK_UPDATE_TEMPLATE'] = CONFIG_BULK_UPDATE_TEMPLATE
    app.jinja_env.globals['MSNOISE_INDEX_TEMPLATE'] = MSNOISE_INDEX_TEMPLATE


    return admin


# ============================================================================
# Main Entry Point
# ============================================================================

def main(port=5099):
    """Main entry point for the MSNoise admin interface"""

    logger.info("Starting MSNoise Admin Interface")

    # Create admin interface
    admin = create_admin_app()

    # Configure Flask app
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

    # Add template filters
    @app.template_filter('datetime')
    def datetime_filter(value):
        if value is None:
            return ''
        return value.strftime('%Y-%m-%d %H:%M:%S')

    # Error handlers
    @app.errorhandler(404)
    def not_found(error):
        return render_template_string('<h1>Page not found</h1>'), 404

    @app.errorhandler(500)
    def internal_error(error):
        db_session.rollback()
        return render_template_string('<h1>Internal server error</h1>'), 500

    # Cleanup
    @app.teardown_appcontext
    def shutdown_session(exception=None):
        db_session.close()

    # Start the application
    try:
        app.run(host='0.0.0.0', port=port, debug=True)
    except KeyboardInterrupt:
        logger.info("Shutting down MSNoise Admin Interface")
    finally:
        db_session.close()


if __name__ == '__main__':
    main()
